from setuptools import setup


setup(
    name='csst_msc_instrument_effect',
    author='xiezhou',
    version='0.1',
    packages=['csst_msc_instrument_effect'],
    entry_points={
        'console_scripts': [
            'csst_msc_instrument_effect=csst_msc_instrument_effect.msc_iec:main'
        ]
    }
)
