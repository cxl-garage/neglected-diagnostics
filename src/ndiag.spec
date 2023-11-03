
# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files
from PyInstaller.utils.hooks import copy_metadata

datas = [("./neglected-diagnostics-venv/Lib/site-packages/streamlit/runtime", "./streamlit/runtime")]
datas += collect_data_files("streamlit")
datas += copy_metadata("streamlit")
datas += [("./app", "./app")]
datas += [("./genetic_testing", "./genetic_testing")]
datas += [("./utils", "./utils")]
datas += [("./neglected-diagnostics-venv/Lib/site-packages/Bio", "./Bio")]
datas += [("./neglected-diagnostics-venv/Lib/site-packages/st_aggrid", "./st_aggrid")]
datas += [("./neglected-diagnostics-venv/Lib/site-packages/pyarrow", "./pyarrow")]

block_cipher = None

a = Analysis(
    ['ndiag.py'],
    pathex=["."],
    binaries=[],
    datas=datas,
    hiddenimports=['Levenshtein', 'fuzzysearch', 'st_aggrid', 'decouple'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='ndiag',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['assets\\CXL.ico'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='ndiag',
)

