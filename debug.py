# with open() as f_mixed, open('log_ax.txt') as f_ax:
for i, f in enumerate(('log_mixed.txt', 'log_ax.txt')):
    with open(f[:-4] + '_norm.txt', 'w') as f_out, open(f) as f_in:
        time = -1
        num = 0
        for line in f_in:
            _, t = line.split(' - ')
            if time == t:
                num += 1
            else:
                if num != 0:
                    f_out.write(f'{int(time) - 25*0} - {num}\n')
                time = t
                num = 1
        f_out.write(f'{int(time) - 25} - {num}\n')
