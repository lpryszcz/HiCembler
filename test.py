def truncate(t):
    for i,n in enumerate(t.traverse(),1):
        if t.get_distance(n, topology_only=1)>3 and len(get_chromosome(n))==1:
            n.name="%s %s leaves"%(get_chromosome(n).pop(), len(n))
            for _n in n.get_children():
                n.remove_child(_n)
    return t
--

