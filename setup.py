"""
Created on 14/8/2014

@author: victor
"""

if __name__ == '__main__': # Compatibility with sphynx
    from setuptools import setup

    setup(
          name='HIVProteaseMutator',
          version='0.3.1',
          author='Victor Alejandro Gil Sepulveda',
          author_email='victor.gil.sepulveda@gmail.com',
          packages=[
                    'hivprotmut',
                    'hivprotmut.sequences',
                    'hivprotmut.mutation',
                    'hivprotmut.external',
                    'hivprotmut.external.blast',
                    'hivprotmut.external.plop',
                    'hivprotmut.external.proteinwizard',
          ],

          install_requires=[
            "ProDy>=1.4.2"
          ]
    )
