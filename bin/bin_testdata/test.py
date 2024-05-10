import glob

def get_scope_bc(chemistry, assets_dir):
    """Return (linker file path, whitelist file path)"""
    try:
        linker_f = glob.glob(f'{assets_dir}/whitelist/{chemistry}/linker*')[0]
        whitelist_f = f'{assets_dir}/whitelist/{chemistry}/bclist'
    except IndexError:
        return None, None
    return linker_f, whitelist_f

a = get_scope_bc('scopeV2.2.1', '/SGRNJ06/randd/USER/liuzihao/work/scsnp/assets')

print(a)