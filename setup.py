from setuptools import setup, find_packages

setup(
    name="bad",
    version="0.0.1",
    packages=find_packages(),
    entry_points={
        'console_scripts': ['bad = app.__main__:main']
    },
    author='Julius Wiegert',
    author_email='julius-wiegert@web.de',
    description='Tool for the estimation of the difficulty of phylogenetic placements',
    long_description='BAD provides a way to estimate the difficulty of phylogenetic placements. It can help to understand why specific placements are easy or difficult due to the usage of Shapley values.',
    include_package_data=True,
    install_requires=[
        "pandas",
        "numpy",
        "ete3",
        "biopython",
        "networkx",
        "scipy",
        "lightgbm",
        "shap",
        "rich-argparse",
        "pyprobables"

    ],
    package_data={
        'bad': ['models/*.pkl'],
    },
)
