# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='rdsolver',
    version='0.0.1',
    author='Shiling Liang',
    author_email='liangshiling42@gmail.com',
    description='simulate reaction diffusion system',
    long_description=open('README.md').read(),
    license='Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)',
    install_requires=['numpy','matplotlib','scipy'],
    url='https://github.com/Shiling42/reaction-diffusion-solver',
    packages=find_packages(exclude=('tests', 'docs')),
    test_suite='tests'
)