"""
OCT文件读取模块
读取自定义.oct文件格式
"""

import numpy as np
import struct
from typing import Tuple, Dict, Any
from pathlib import Path


class OCTFileReader:
    """OCT文件读取器"""
    
    def __init__(self, filename: str):
        """
        初始化OCT文件读取器
        
        Args:
            filename: OCT文件路径
        """
        self.filename = filename
        self.fid = None
        self.header = {}
        
    def open(self):
        """打开文件"""
        self.fid = open(self.filename, 'rb')
        
    def close(self):
        """关闭文件"""
        if self.fid:
            self.fid.close()
            self.fid = None
    
    def __enter__(self):
        self.open()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
    
    def read_header(self) -> Dict[str, Any]:
        """
        读取文件头信息
        
        Returns:
            包含文件头信息的字典
        """
        if not self.fid:
            self.open()
        
        self.fid.seek(0)
        
        # 读取基本参数
        self.header['bob'] = struct.unpack('I', self.fid.read(4))[0]  # uint32
        self.header['SPL'] = struct.unpack('d', self.fid.read(8))[0]  # double
        self.header['nX'] = struct.unpack('I', self.fid.read(4))[0]   # uint32, A-lines数量
        self.header['nY'] = struct.unpack('I', self.fid.read(4))[0]   # uint32, B-scans数量
        self.header['Boffset'] = struct.unpack('I', self.fid.read(4))[0]
        self.header['Blength'] = struct.unpack('I', self.fid.read(4))[0] + 1
        self.header['Xcenter'] = struct.unpack('d', self.fid.read(8))[0]
        self.header['Xspan'] = struct.unpack('d', self.fid.read(8))[0]
        self.header['Ycenter'] = struct.unpack('d', self.fid.read(8))[0]
        self.header['Yspan'] = struct.unpack('d', self.fid.read(8))[0]
        self.header['frame_per_pos'] = struct.unpack('I', self.fid.read(4))[0]  # B-scan重复次数
        self.header['n_dataset'] = struct.unpack('I', self.fid.read(4))[0]      # 体积扫描重复次数
        self.header['ProtMode'] = struct.unpack('I', self.fid.read(4))[0]
        
        # 跳过4字节 (v10)
        self.fid.seek(4, 1)
        
        # 读取背景数据1
        sizeBck1 = struct.unpack('I', self.fid.read(4))[0]
        self.header['Bck1'] = np.frombuffer(self.fid.read(sizeBck1 * 2), dtype=np.int16)
        
        # 读取KES1
        sizeKES1 = struct.unpack('II', self.fid.read(8))
        KES1_data = np.frombuffer(self.fid.read(sizeKES1[1] * 8), dtype=np.float64)
        self.header['KES1'] = KES1_data * sizeKES1[1]
        
        # 读取背景数据2
        sizeBck2 = struct.unpack('I', self.fid.read(4))[0]
        self.header['Bck2'] = np.frombuffer(self.fid.read(sizeBck2 * 2), dtype=np.int16)
        
        # 读取KES2
        sizeKES2 = struct.unpack('II', self.fid.read(8))
        KES2_data = np.frombuffer(self.fid.read(sizeKES2[1] * 8), dtype=np.float64)
        self.header['KES2'] = KES2_data * sizeKES2[1]
        
        # 读取色散系数
        self.header['disp_coef'] = struct.unpack('d', self.fid.read(8))[0]
        
        # 计算衍生参数
        self.header['nR'] = self.header['frame_per_pos']
        self.header['IMGheight'] = self.header['Blength'] // 2
        self.header['nY_actual'] = self.header['nY'] // self.header['nR']
        
        return self.header
    
    def read_bscan(self, iY: int, nr: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        读取单个B-scan的数据
        
        Args:
            iY: B-scan索引 (0-based)
            nr: 重复次数
        
        Returns:
            Bs1, Bs2: 两个通道的数据 [Blength×nX×nr]
        """
        if not self.fid:
            self.open()
        
        bob = self.header['bob']
        SPL = int(self.header['SPL'])
        nX = self.header['nX']
        nR = self.header['nR']
        Blength = self.header['Blength']
        
        Bs1 = np.zeros((Blength, nX, nr), dtype=np.float64)
        Bs2 = np.zeros((Blength, nX, nr), dtype=np.float64)
        
        # 计算文件偏移
        offset = bob + (SPL * nX * 2 + 2) * 2 * (nR * iY)
        self.fid.seek(offset)
        
        for ic in range(nr):
            # 跳过4字节
            self.fid.seek(4, 1)
            
            # 读取数据
            B = np.frombuffer(self.fid.read(SPL * nX * 2 * 2), dtype=np.int16)
            B = B.reshape((SPL, nX * 2), order='F').astype(np.float64)
            
            Bs1[:, :, ic] = B[:Blength, :nX]
            Bs2[:, :, ic] = B[:Blength, nX:]
        
        return Bs1, Bs2
    
    def read_reference(self, n_ref: int = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        读取参考数据
        
        Args:
            n_ref: 参考帧数，默认自动计算
        
        Returns:
            Ref_ch1, Ref_ch2: 参考数据
        """
        if not self.fid:
            self.open()
        
        bob = self.header['bob']
        SPL = int(self.header['SPL'])
        nX = self.header['nX']
        nY = self.header['nY_actual']
        
        if n_ref is None:
            n_ref = min(50000, (nY * nX) // 2 * 2) // nX
        
        self.fid.seek(bob)
        
        Ref_ch1 = np.zeros((SPL, nX, n_ref), dtype=np.float64)
        Ref_ch2 = np.zeros((SPL, nX, n_ref), dtype=np.float64)
        
        for i_ref in range(n_ref):
            self.fid.seek(4, 1)
            BB = np.frombuffer(self.fid.read(SPL * nX * 2 * 2), dtype=np.int16)
            BB = BB.reshape((SPL, nX * 2), order='F').astype(np.float64)
            Ref_ch1[:, :, i_ref] = BB[:, :nX]
            Ref_ch2[:, :, i_ref] = BB[:, nX:]
        
        # 计算平均参考
        Ref_ch1 = np.tile(np.mean(np.mean(Ref_ch1, axis=1, keepdims=True), axis=2, keepdims=True), (1, nX, 1))[:, :, 0]
        Ref_ch2 = np.tile(np.mean(np.mean(Ref_ch2, axis=1, keepdims=True), axis=2, keepdims=True), (1, nX, 1))[:, :, 0]
        
        return Ref_ch1, Ref_ch2


def load_oct_file(filename: str) -> Tuple[Dict[str, Any], OCTFileReader]:
    """
    便捷函数：加载OCT文件并返回头信息和读取器
    
    Args:
        filename: OCT文件路径
    
    Returns:
        header: 文件头信息
        reader: OCT文件读取器对象
    """
    reader = OCTFileReader(filename)
    header = reader.read_header()
    return header, reader


if __name__ == "__main__":
    import sys
    
    # 测试用例
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        
        print(f"读取文件: {filename}")
        
        with OCTFileReader(filename) as reader:
            header = reader.read_header()
            
            print("\n文件头信息:")
            print(f"  SPL: {header['SPL']}")
            print(f"  nX: {header['nX']}")
            print(f"  nY: {header['nY']}")
            print(f"  nY_actual: {header['nY_actual']}")
            print(f"  Blength: {header['Blength']}")
            print(f"  frame_per_pos: {header['frame_per_pos']}")
            print(f"  disp_coef: {header['disp_coef']}")
    else:
        print("用法: python oct_file_reader.py <oct文件路径>")
