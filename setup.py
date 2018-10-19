from setuptools import setup

setup(
        name='dss_workflow',
        version='0.1',
        author='Stephen Kraemer',
        author_email='stephenkraemer@gmail.com',
        license='MIT',
        packages=['dss_workflow'],
        install_requires=[
            'pandas',
        ],
        python_requires='>=3.6',
)
