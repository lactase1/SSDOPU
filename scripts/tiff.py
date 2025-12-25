import os
import shutil

src_root = r'D:/1-Liu Jian/yongxin.wang/PSOCT/2025-9-19'
dst_root = os.path.join(src_root, 'TIFF')
os.makedirs(dst_root, exist_ok=True)

for param_folder in os.listdir(src_root):
    param_path = os.path.join(src_root, param_folder)
    if not os.path.isdir(param_path):
        continue
    for data_folder in os.listdir(param_path):
        data_path = os.path.join(param_path, data_folder)
        tiff_path = os.path.join(data_path, 'tiff')
        if os.path.isdir(tiff_path):
            for fname in os.listdir(tiff_path):
                if fname.lower().endswith('.tiff'):
                    # 类型名为最后一个下划线后所有内容
                    parts = fname.split('_')
                    if len(parts) < 2:
                        continue
                    # 类型名 = 最后一个下划线后所有内容
                    type_name = '_'.join(parts[-3:])  # 适配你这种命名，通常为3段
                    # 也可以用 type_name = '_'.join(parts[2:]) 取第3段及以后
                    # 但更保险的做法是：类型名 = 原文件名去掉前缀（即第一个参数目录和数据目录相关的部分）
                    # 你实际想要的是：找到第一个出现“_1-”或“_2-”等的地方
                    for i, p in enumerate(parts):
                        if p.startswith(('1-', '2-', '3-', '4-', '5-', '6-', '7-', '8-', '9-')):
                            type_name = '_'.join(parts[i:])
                            break
                    else:
                        type_name = fname  # fallback

                    # 目标子文件夹
                    type_dir = os.path.join(dst_root, type_name.replace('.tiff',''))
                    os.makedirs(type_dir, exist_ok=True)
                    # 新文件名
                    new_name = f"{param_folder}_{type_name}"
                    src_file = os.path.join(tiff_path, fname)
                    dst_file = os.path.join(type_dir, new_name)
                    shutil.copy2(src_file, dst_file)
                    print(f"复制: {src_file} -> {dst_file}")

print("全部完成！")