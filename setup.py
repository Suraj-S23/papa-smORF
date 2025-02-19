from setuptools import setup, find_packages

setup(
    name="papa_smorf",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'pandas>=2.0.0',
        'numpy>=1.24.0',
        'pyyaml>=6.0.1',
        'biopython>=1.81',
        'matplotlib>=3.5.0',
        'seaborn>=0.11.0',
        'tqdm>=4.65.0',
        'pytest>=7.0.0',
        'psutil>=5.9.0'  
    ]
)
