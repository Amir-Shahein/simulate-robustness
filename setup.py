from setuptools import setup
import robustness_analysis

setup(
	name='robustness_analysis',
    version=robustness_analysis.__version__,
    description='robustness-analysis',
    url='none',
    author='AmirShahein',
    author_email='amir.shahein@epfl.ch',
    license='none',
    packages=['robustness_analysis'],
    python_requires='>3.6',
    install_requires=['numpy', 'pandas'],
    zip_safe=False
)
