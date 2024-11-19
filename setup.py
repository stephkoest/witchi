from setuptools import setup, find_packages

setup(
    name='witchi',
    version='0.1.0-alpha',
    description='A compositional bias pruning tool for multiple sequence alignments.',
    author='Stephan Koestlbacher',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'witchi=witchi.witchi:main',
        ],
    },
    install_requires=[
        'numpy',
        'scipy',
        'biopython',
        'joblib'
    ],
)
