from setuptools import setup, find_packages

setup(
    name='mgwas',
    version='0.1.0',
    author='Zeyang Shen',
    author_email='zeyang.shen@wsu.edu',
    description='A microbial GWAS tool for trait association from shotgun metagenomic data',
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    url='https://github.com/skinmicrobiome/VITALITY_skin_atopy',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[],
    entry_points={
        'console_scripts': [
            'mgwas=mgwas.cli:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.12',
)