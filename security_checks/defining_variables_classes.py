class Prueba:

    _value1 = 0

    def __init__(self, value1=_value1):
        self.value1 = value1

    def reset(self):
        self.value1 = self._value1

a = Prueba()

print(a.value1)

a.value1 = 3

a.reset()

print(a.value1)

a.value1 = 3

print(a.value1)
