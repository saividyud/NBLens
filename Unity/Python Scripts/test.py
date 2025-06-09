print('Hello world')

import platform
import psutil

print(platform.machine())
print(platform.version())
print(platform.system())
print(platform.processor())

memory = psutil.virtual_memory()

print(memory.total / (1024 ** 3))
print(memory.used / (1024 ** 3))
print(memory.available / (1024 ** 3))
