obj, a = 'gp123',10.123413549531298956
print(f'{a:.2f}')
print(f'{a:.3e}')
print('{:.8f}'.format(a))
print('{:.3e}'.format(a))
combined_title = f'Galaxy {obj} STARLIGHT synthesis mass fraction'\
                + '\n' \
                + r'$Log(M_{{\star}})={:.2f}$'.format(a)#
print(combined_title)
b = ('This is my '
      'Name')
print(b, type(b))
c = b + ('Vital the god')
print(c)

