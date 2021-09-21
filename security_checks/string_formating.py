obj, a = 'gp123', 10.583413549531298956
print(f'{a:.2f}')
print(f'{a:.3e}')
print(f'{a:.0f}')
print('{:.8f}'.format(a))
print('{:.3e}'.format(a))
print('{:.4e}'.format(a))

cifras = 4
print(f'{a:.{cifras}e}')
print(f'{a:.{cifras}f}')
print(r'$Log(M_{{\star}})=\num{{{:.3e}}}$'.format(a))


combined_title = f'Galaxy {obj} STARLIGHT synthesis mass fraction'\
                + '\n' \
                + r'$Log(M_{{\star}})={:.2f}$'.format(a)#
print(combined_title)
b = ('This is my '
      'Name')
print(b, type(b))
c = b + ('Vital the god')
print(c)

a = 'H1_6563A_w1_center'
print(a[a.rfind('_')+1:])

