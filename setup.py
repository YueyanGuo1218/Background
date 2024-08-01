from setuptools import setup, find_packages

setup(
    name='Background',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pysam',
        'numpy'
    ],
    entry_points={
        'console_scripts': [
            'ssbg=ssbg.main:main',
        ],
    },
    author='Yueyang Guo',
    author_email='pot_fe@outlook.com',
    description='A command-line tool for calculate background for Strand-Seq BAM files',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/YueyanGuo1218/Background',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.6',
)
