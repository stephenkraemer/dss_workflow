from setuptools import setup

setup(
        name='dmr_calling',
        version='0.1',
        author='Stephen Kraemer',
        author_email='stephenkraemer@gmail.com',
        license='MIT',
        packages=['dmr_calling'],
        install_requires=[
            'pandas',
        ],
        python_requires='>=3.6',
)
