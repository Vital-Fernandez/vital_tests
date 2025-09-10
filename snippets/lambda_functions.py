a = f"{(lambda x: x * 37) (2)}"
print(a)
# map(lambda x: 'L({})'.format(x[x.find('_') + 1:len(x) - 1]), line_comps)
b = [1000, 2000, 3000]
ref = "+".join(map("L({:.0f})".format, b))
print(ref)
ion = 'He1'
ion = ion if ion not in ('He1', 'He2') else ion + 'r'
print(ion)