class Car(object):
    condition = "new"

    def __init__(self, model='a', color='b', mpg='c'):
        self.model = model
        self.color = color
        self.mpg = mpg

    def __str__(self):
        return '{},{},{}'.format(self.model, self.color, self.mpg)

class ElectricCar(Car):

    def __init__(self, battery_type, model, color, mpg):
        self.battery_type=battery_type
        super().__init__(model, color, mpg)

    def __str__(self):
        return '{}, {},{},{}'.format(self.battery_type, self.model, self.color, self.mpg)

    def reset(self):
        super().__init__()


bmv1 = ElectricCar('gm', 'bmv', 'red', 'coso')
print(bmv1)
bmv1.reset()
print(bmv1)

