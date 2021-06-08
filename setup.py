# -*- coding: utf-8 -*-
"""\
Allison MacLeay

"""

import re
from distutils.core import setup
from setuptools import find_packages


def read_requirements():
    """Reads the list of required packages from requirements.txt"""
    def parse_one(line):
        """Parses a single requirement"""
        line = line.strip()
        m = re.match(r'-e.*/(.*)@.*#egg=.*', line)  # for git links: use the package name
        if m is not None:
            return m.group(1)
        else:
            return line

    with open('requirements.txt', 'r') as f:
        return [parse_one(line) for line in f if '://' not in line]


setup(
    name='biollama',
    version='0.0.2',
    packages=find_packages(),
    install_requires=read_requirements(),
    url='',
    license='MIT',
    author='Allison MacLeay',
    author_email='allison.macleay@gmail.com',
    description='bioinformatics tools'
)
