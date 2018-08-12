#Testing For Loop Else

tmp = [1,2,3,4]
for i in tmp:
    print i
    continue
else:
    print 'ELSE1'

tmp = [1,2,3,4]
for i in tmp:
    print i
    if i == 2:
        break
else:
    print 'ELSE2'

