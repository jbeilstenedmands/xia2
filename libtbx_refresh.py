import dials.precommitbx.nagger

dials.precommitbx.nagger.nag()


def _install_xia2_setup():
    """Install xia2 as a regular/editable python package"""
    import subprocess
    import sys

    import libtbx.load_env

    # Call pip
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-deps",
            "-e",
            libtbx.env.dist_path("xia2"),
        ],
        check=True,
    )


def _show_xia2_version():
    try:
        import xia2

        print(xia2.__version_string__)
    except ModuleNotFoundError:
        print("Can't tell xia2 version")


_install_xia2_setup()
_show_xia2_version()
