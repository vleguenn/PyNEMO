Troubleshooting
===============

1. pyNEMO crashing in MacOSX (Yosemite)?

*  Downgrade the scipy package to 0.15

2. How to make pyNEMO to work behind firewall/proxy?

*  Set the environment variable http_proxy. eg. in Linux export http_proxy=<proxy-server>:<proxy-port>

3. Getting this error 'Warning: Please make sure pyjnius is installed and jvm.dll/libjvm.so/libjvm.dylib is in the path' ?

*  This error is displayed when the application cannot find the java installation on the local machine. please install a java 7.x runtime from http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html and append the path to the library in the system path. eg. on windows SET PATH="C:\\Program Files (x86)\\Java\\jre1.7\\bin\\client"
