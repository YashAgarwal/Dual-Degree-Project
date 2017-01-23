import glob,os.path
filesDepth3 = glob.glob('*/*/*')
dirs = filter(lambda f: os.path.isdir(f), filesDepth3)
#print dirs
variables = ['c']
for var in variables:
    print "Printing variable %s" % var
    for path in dirs:
        f = open("%s/README" % (path))
        a = f.read()
        a = a.split('\n')
        #N1 = a[a.find("time_steps") + len("time_steps") + 3]
        #print N1
        for b in a:
            if "time_steps" in b and "M0" not in b:
                #print b
                c = b.split(' = ')
                N1 = int(c[1])
            elif "time_steps_M0" in b:
                #print b
                c = b.split(' = ')
                N2 = int(c[1])
            elif "file_timer" in b:
                c = b.split(' = ')
                S = int(c[1])
            elif "n_x" in b:
                c = b.split(' = ')
                N = int(c[1])
            else:
                #print b
                pass

        print "N1 = %d" % N1
        print "N2 = %d" % N2
        print "S = %d" % S
        print "N = %d" % N

        print path
        os.system("mkdir -p %s/analysis/%s_profile" % (path, var))
        cmd = "echo \"outputPath='%s'; name='%s'; N1=%d; N2=%d; S=%d;N=%d; load 'makeMovie.gnu';\" | gnuplot\n" % (path, var, N1, N2, S, N)
        os.system(cmd)
print "Done"
