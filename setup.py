# read the contents of your README file
from pathlib import Path

from setuptools import find_packages, setup

print("I'm here")
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="raytrace",
    version="0.1",
    packages=find_packages(),
    package_dir={"rt": "rt"},
    package_data={"rt": []},
    data_files=[("rt", [])],
    include_package_data=True,
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    author="Shibaji Chakraborty",
    author_email="chakras4@erau.edu",
    maintainer="Shibaji Chakraborty",
    maintainer_email="chakras4@erau.edu",
    license="GNU GPL License",
    description=long_description,
    long_description=long_description,
    install_requires=[],
    keywords=["python", "raytracing using Pharlap"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    url="https://github.com/shibaji7/raytrace",
)
