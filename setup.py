from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="bad-phylo",
    version="0.1.0",
    packages=find_packages(),
    entry_points={
        'console_scripts': ['bad = app.__main__:main']
    },
    author='Julius Wiegert',
    author_email='julius-wiegert@web.de',
    description='Tool for the estimation of the difficulty of phylogenetic placements',
    long_description=long_description,
    long_description_content_type='text/markdown',
    include_package_data=True,
    install_requires=[
        "pandas",
        "numpy",
        "ete3",
        "biopython",
        "networkx",
        "scipy",
        "lightgbm==4.1.0",
        "shap",
        "rich-argparse",
        "pyprobables",
        "matplotlib"

    ],
    package_data={
        'bad': ['models/*.pkl'],
    },
)
