from setuptools import setup

setup(
        name='dss_workflow',
        version='0.1',
        author='Stephen Kraemer',
        author_email='stephenkraemer@gmail.com',
        license='MIT',
        packages=['dss_workflow'],
        package_data={
            '': ['*.R', '*.snakefile', '*.yml', '*.sh'],
        },
        install_requires=[
            'pandas>=0.23',
        ],
        python_requires='>=3.6',
)
