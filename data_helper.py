import os
import shutil

# 原始目录和目标目录
src_dir = r"D:\1-Liu Jian\yongxin.wang\Data"
dst_dir = r"D:\1-Liu Jian\yongxin.wang\SamplesData"

os.makedirs(dst_dir, exist_ok=True)

# 检查源目录是否存在
print(f"源目录: {src_dir}")
print(f"目录是否存在: {os.path.exists(src_dir)}")
print(f"目标目录: {dst_dir}")

# 用字典保存每种类型的第一个文件
selected = {}

oct_files_found = 0
for root, _, files in os.walk(src_dir):
    for f in files:
        if f.endswith(".oct"):
            oct_files_found += 1
            # 去掉时间戳：从第一个 '.' 后面开始取，或者从第一个 '_' 后面开始取
            # 例如: 2024.03.29_19.10.06_1.blueBarDisk_QWP_180_sample_90.oct
            # 去掉前面的时间戳，得到 blueBarDisk_QWP_180_sample_90.oct
            parts = f.split("_", 2)
            key = parts[-1] if len(parts) >= 3 else f  # 如果不满足就保留原名

            # 只保留第一个出现的同类文件
            if key not in selected:
                selected[key] = os.path.join(root, f)

print(f"遍历过程中找到 {oct_files_found} 个 .oct 文件")
print(f"找到 {len(selected)} 个不同类型的 OCT 文件，准备复制...")

for key, fullpath in selected.items():
    dst_path = os.path.join(dst_dir, os.path.basename(fullpath))
    if os.path.exists(dst_path):
        print(f"跳过已有文件: {dst_path}")
    else:
        print(f"复制: {fullpath} -> {dst_path}")
        shutil.copy2(fullpath, dst_path)

print("全部处理完成")
