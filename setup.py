from setuptools import setup, find_packages

setup(
    name = 'fempy',
    version = '1.0.0',
    description = 'Yet another FEM suite for 2D structural analysis',
    author = 'Nick An',
    author_email = 'nickanpers@gmail.com',
    url = 'https://github.com/nikolisan/FEMpython/',
    install_requires = [
        "numpy",
        "matplotlib",
        "scipy",
        "triangle",
        "loguru"
    ],
    packages=find_packages()
)