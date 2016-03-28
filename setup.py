from setuptools import setup, find_packages
def calculate_version(inputfile):
    version_list =  [x.split('\'')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

package_version = calculate_version('./pyutilities/pyutilities.py')

setup(
    name='PyUtilities',
    version=package_version,
    #py_modules=['pyutilities'],
    pacages=find_packages(),
    install_requires=[
        'Click', 'beautifulsoup4', 'biopython',
    ],
    entry_points="""
    [console_scripts]
    pyutilities=pyutilities.pyutilities:cli
    """
)
