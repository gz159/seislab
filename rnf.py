import os
import glob

# 获取当前目录下所有 .mseed 文件
files = glob.glob('*.mseed')

for old_name in files:
    # 去掉扩展名，分割文件名
    name_without_ext = os.path.splitext(old_name)[0]
    parts = name_without_ext.split('_')
    
    # 检查是否有足够的组成部分
    if len(parts) >= 8:
        # 获取第5、6、8部分 (索引是 4, 5, 7)
        new_base_name = f"{parts[4]}_{parts[5]}_{parts[7]}"
        new_name = new_base_name + '.mseed'
        
        # 重命名文件
        os.rename(old_name, new_name)
        print(f"重命名: {old_name} -> {new_name}")
    else:
        print(f"跳过 '{old_name}'：文件名部分不足8个。")