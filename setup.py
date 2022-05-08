# Setup File
import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
# "packages": ["os"] is used as example only
build_exe_options = {"packages": ["os"], "include_files": [ "images\Aircraft_image.png", "images\Ramjet.png",
                                                            "images\Turbojet_w_AB.png", "images\Turbojet_wo_AB.png",
                                                            "images\Turbofan_w_AB.png", "images\Turbofan_wo_AB.png",
                                                            "images\Mixed_flow_w_AB.png", "images\Mixed_flow_wo_AB.png",
                                                           "dependencies\Turbofan.py", "dependencies\Mixed_flow_Turbofan.py"],
                     "excludes": [""]}

# base="Win32GUI" should be used only for Windows GUI app
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(
    name="Aircraft Engine Design",
    version="1.1",
    description="Aircraft Engine Design",
    options={"build_exe": build_exe_options},
    executables=[Executable("main.py", base=base)],
)