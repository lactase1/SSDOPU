#!/usr/bin/env python3
"""
reduce_images.py

按指定步长（默认 5）从每个子目录的图片文件中保留第 1、6、11... 张图片。
支持模式：copy（默认，将保留的图片复制到子目录下的 `reduced` 文件夹），
move（将保留的图片移动到 `reduced`），
inplace（删除不保留的图片，需要 --yes 确认）。

示例：
python reduce_images.py -r "C:\yongxin.wang\Data\Enface_struct\2024.09.04_18.42.53_1.SLnoBG_disk\struct_pngs" --mode copy --dry-run
python reduce_images.py -r "C:\yongxin.wang\Data\Enface_struct\2024.09.04_18.42.53_1.SLnoBG_disk\struct_pngs" --mode inplace --yes

"""

from pathlib import Path
import argparse
import shutil
import re
import sys

IMAGE_EXTS = {'.png', '.jpg', '.jpeg', '.bmp', '.tif', '.tiff'}

def natural_key(s: str):
    # 分割字符串中的数字，使排序更符合人类习惯
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]


def process_folder(folder: Path, args):
    files = [p for p in sorted(folder.iterdir(), key=lambda p: natural_key(p.name)) if p.is_file() and p.suffix.lower() in args.exts]
    n = len(files)
    if n == 0:
        return {'processed': False, 'reason': 'no_images', 'n': 0}
    if n < args.keep_step:
        return {'processed': False, 'reason': 'too_few', 'n': n}

    kept_idx = [i for i in range(n) if (i - args.start_index) % args.keep_step == 0]
    kept = [files[i] for i in kept_idx]
    removed = [files[i] for i in range(n) if i not in kept_idx]

    if args.mode in ('copy', 'move'):
        dest = folder / args.dest_name
        if not args.dry_run:
            dest.mkdir(parents=True, exist_ok=True)
        for f in kept:
            destf = dest / f.name
            if args.dry_run:
                print(f"[DRY-RUN] {args.mode} {f} -> {destf}")
            else:
                if args.mode == 'copy':
                    shutil.copy2(f, destf)
                else:
                    shutil.move(str(f), str(destf))
        return {'processed': True, 'n': n, 'kept': len(kept), 'removed': 0 if args.mode != 'inplace' else len(removed)}

    if args.mode == 'inplace':
        if args.dry_run:
            for f in removed:
                print(f"[DRY-RUN] delete {f}")
            return {'processed': True, 'n': n, 'kept': len(kept), 'removed': len(removed), 'dry': True}
        if not args.yes:
            return {'processed': False, 'reason': 'confirm_required', 'n': n}
        # 执行删除
        for f in removed:
            try:
                f.unlink()
            except Exception as e:
                print(f"Failed to delete {f}: {e}")
        return {'processed': True, 'n': n, 'kept': len(kept), 'removed': len(removed)}

    return {'processed': False, 'reason': 'unknown_mode'}


def main():
    p = argparse.ArgumentParser(description='按固定步长从每个子文件夹中保留图片（默认每5张保留1张）。')
    p.add_argument('-r', '--root', type=Path, default=Path('.'), help='顶层目录（递归查找包含图片的子目录）')
    p.add_argument('--keep-step', type=int, dest='keep_step', default=5, help='每隔多少张保留一张（默认 5）')
    p.add_argument('--start-index', type=int, default=0, help='起始索引（0 起）以决定保留哪张，默认 0 即保留第1张')
    p.add_argument('--mode', choices=['copy', 'move', 'inplace'], default='copy', help='操作模式：copy/move/inplace（默认 copy）')
    p.add_argument('--dest-name', default='reduced', help='当使用 copy/move 时，保留图片放到的子文件夹名（默认 reduced）')
    p.add_argument('--dry-run', action='store_true', help='仅打印将要执行的操作，不做实际改变')
    p.add_argument('--yes', action='store_true', help='当 --mode inplace 时，确认执行删除（谨慎）')
    p.add_argument('--exts', default=','.join(sorted(IMAGE_EXTS)), help='允许的图片扩展名，逗号分隔（例如 .png,.jpg）')

    args = p.parse_args()
    args.exts = set(ext.strip().lower() if ext.strip().startswith('.') else f'.{ext.strip().lower()}' for ext in args.exts.split(','))
    args.start_index = int(args.start_index)
    args.keep_step = int(args.keep_step)

    if args.keep_step <= 0:
        print('keep-step 必须是正整数。')
        sys.exit(1)
    if args.start_index < 0:
        print('start-index 必须是 >= 0 的整数。')
        sys.exit(1)

    root = args.root
    if not root.exists():
        print(f'指定路径不存在: {root}')
        sys.exit(1)

    # 递归查找包含图片的目录并处理
    folders = []
    for d in root.rglob('*'):
        if d.is_dir():
            # 判断此目录是否包含图片文件
            if any((f.suffix.lower() in args.exts) for f in d.iterdir() if f.is_file()):
                folders.append(d)

    if not folders:
        print('未在子目录中发现任何图片文件。')
        sys.exit(0)

    total = {'folders': 0, 'files': 0, 'kept': 0, 'removed': 0}
    for folder in sorted(folders):
        res = process_folder(folder, args)
        if res.get('processed'):
            print(f"Processed folder: {folder} | total {res.get('n')} => kept {res.get('kept',0)} removed {res.get('removed',0)}")
            total['folders'] += 1
            total['files'] += res.get('n',0)
            total['kept'] += res.get('kept',0)
            total['removed'] += res.get('removed',0)
        else:
            reason = res.get('reason')
            if reason == 'no_images':
                pass
            elif reason == 'too_few':
                print(f"Skip folder (fewer than keep-step): {folder} (found {res.get('n')})")
            elif reason == 'confirm_required':
                print(f"Skipped deleting in {folder}: run with --yes to confirm inplace deletion.")
            else:
                print(f"Skipped folder {folder}: {res}")

    print('\nSummary:')
    print(f"Folders processed: {total['folders']}")
    print(f"Total files seen: {total['files']}")
    print(f"Total kept: {total['kept']}")
    print(f"Total removed (or moved out): {total['removed']}")

if __name__ == '__main__':
    main()
