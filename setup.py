from setuptools import setup
def readme():
    with open('README.md') as f:
        README = f.read()
    return README


setup(
    name="geomeshconv",
    version="0.0.9",
    description="A Python package to convert the geology to in mesh generator (GMSH).",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/Ali1990dashti/GeoMeshPy",
    author="Ali Dashti",
    author_email="Ali.dashti@kit.edu",
    license="MIT",
    packages=["GeoMeshPy"],
    include_package_data=True,
    install_requires=["numpy", "pandas", "numpy_indexed"],
)
