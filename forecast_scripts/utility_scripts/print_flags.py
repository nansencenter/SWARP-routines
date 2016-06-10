f     = open('flags','r')
lines = f.readlines()
f.close()

flags = []
for lin in lines:
   lst      = lin.split()
   ignore   = False
   for ss in lst:
      if ss=='#':
         ignore   = True # ignore rest of line (comment)
      elif not ignore:
         if ('-D' in ss) and ('$' not in ss):
            print(ss)
