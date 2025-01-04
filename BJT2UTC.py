from datetime import datetime, timedelta, timezone

# 给定的北京时间字符串
beijing_time_str = "2025-01-04 06:45:57"

# 定义格式
fmt = "%Y-%m-%d %H:%M:%S"

# 创建一个带有时区信息的 datetime 对象 (UTC+8)
beijing_dt = datetime.strptime(beijing_time_str, fmt)

# 添加时区信息（UTC+8）
beijing_tz = timezone(timedelta(hours=8))
beijing_dt = beijing_dt.replace(tzinfo=beijing_tz)

# 转换为 UTC 时间
utc_dt = beijing_dt.astimezone(timezone.utc)

# 打印结果
print("UTC time:", utc_dt.strftime(fmt))