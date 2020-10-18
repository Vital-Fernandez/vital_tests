def operations(function):
    def coso(*args, **kwargs):
        try:
            print("Numbers:", *args)
            response = function(*args, **kwargs)
            print("Result:", response)
            return response
        except Exception as err:
            print(err)

    return coso


@operations
def subtraction(num_a, num_b, bicho=None):
    print(bicho)
    return num_a - num_b


@operations
def division(num_a, num_b):
    return num_a / num_b


subtraction(5, 2, bicho='Coso')
print()
division(5, 0)