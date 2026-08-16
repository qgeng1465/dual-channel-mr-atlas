#!/usr/bin/env python3
# =============================================================================
# M11_find_tar_offset.py — 用 Range 请求顺序跳读 tar 头，定位目标染色体文件偏移
# =============================================================================
# 原理：标准 ustar 每个条目 = 512B 头 + ceil(size/512)*512B 数据。
#   顺序取 512B 头，解析 name/size 算下一偏移，直到找到目标条目。
#   全程只下载每个条目的 512B 头（~30 个小请求），不下载数据。
# 用法：python3 M11_find_tar_offset.py <synid> <目标chr字符串>
#   e.g. python3 M11_find_tar_offset.py syn51468818 chr1
# 输出：目标条目 {offset} {size} {name}
# =============================================================================
import json, os, re, sys, urllib.request as u

synid = sys.argv[1]
wanted = sys.argv[2]
pat = re.compile(rf'_chr{wanted}_')  # 精确染色体（chr1 不匹配 chr16）
tok = os.popen("awk -F' = ' '/authToken/{print $2}' ~/.synapseConfig").read().strip()
proxy = u.ProxyHandler({'https': 'http://127.0.0.1:7890'})
opener = u.build_opener(proxy)

def get_range(offset, length, url):
    req = u.Request(url, headers={'Range': f'bytes={offset}-{offset+length-1}'})
    return opener.open(req, timeout=25).read()

# 取 pre-signed URL
FH = json.load(opener.open(u.Request(
    f'https://repo-prod.prod.sagebase.org/repo/v1/entity/{synid}',
    headers={'Authorization': 'Bearer ' + tok}), timeout=25)).get('dataFileHandleId')
batch = {"includeFileHandles": True, "includePreSignedURLs": True,
         "requestedFiles": [{"fileHandleId": FH, "associateObjectId": synid,
                             "associateObjectType": "FileEntity"}]}
url = json.load(opener.open(u.Request(
    'https://file-prod.prod.sagebase.org/file/v1/fileHandle/batch',
    data=json.dumps(batch).encode(), headers={'Authorization': 'Bearer ' + tok,
                                              'Content-Type': 'application/json'}), timeout=25))['requestedFiles'][0]['preSignedURL']

def read_name(h):
    return h[:100].rstrip(b'\x00').decode('utf-8', 'replace')

def read_size(h):
    field = h[124:136]
    # GNU base-256 编码（首字节最高位为 1，用于 >8GB 条目）
    if field and field[0] & 0x80:
        v = field[0] & 0x7F
        for c in field[1:]:
            v = (v << 8) | c
        return v
    # 标准 ustar：size 字段为八进制（前导 0 + NUL/空格补齐），必须按 8 进制解析
    try:
        return int(field.strip(b'\x00 ').decode('ascii') or '0', 8)
    except (ValueError, UnicodeDecodeError):
        return 0

off = 0
found = None
entries = []
for step in range(200):
    h = get_range(off, 512, url)
    if len(h) < 512:
        break
    name = read_name(h)
    if not name:
        break
    sz = read_size(h)
    typeflag = h[156:157]
    entries.append((off, name, sz, typeflag))
    if typeflag != b'5' and typeflag != b'0' and typeflag != b'':
        pass
    if pat.search(name) and typeflag not in (b'5',):
        found = (off, sz, name)
        break
    # 跳过数据块到下一头
    data_blocks = ((sz + 511) // 512) * 512
    off += 512 + data_blocks

print(f"# {synid} entries_scanned={len(entries)}")
if found:
    print(f"{found[0]} {found[1]} {found[2]}")
else:
    print("# NOT_FOUND. 目录前 8 条目:")
    for off, name, sz, tf in entries[:8]:
        print(f"#   off={off} size={sz} type={tf!r} {name}")
