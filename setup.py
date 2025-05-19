from setuptools import setup

setup(
    name='oad',
    version='1.0',
    description='Simulate oad model',
    author='Tomas Scagliarini',
    author_email='tomscag92@gmail.com',
    packages=['oad'], 
    python_requires=">=3.12",
    install_requires=[
        "numpy",
        "networkx>=3.4",
    ],
)
