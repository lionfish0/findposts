import numpy as np

def find_post_code(img,markerloc,angres=120,sizeres=3):
    optimum_result = None
    for ang in np.linspace(0,np.pi*2,angres):
        dists = np.linspace(0,300,301)
        try:
            vs = img[(np.cos(ang)*dists+markerloc[0]).astype(int),(np.sin(ang)*dists+markerloc[1]).astype(int)]
        except IndexError:
            continue #this angle goes outside image...
        maxval = np.max(vs[10:])
        mid = np.max(vs[10:])/2+np.min(vs[10:])/2
        
        vs[vs>maxval] = maxval #trying to remove retroreflector spike!
        #print(mid,vs)
        pat = np.repeat(np.array([255,0,255,0,0,255])[:,None],20,1).flatten()
        #pat = np.repeat(np.array([255,255,255,255,255,255])[:,None],20,1).flatten()
        norm_vs = vs - np.mean(vs)
        
        maxv, max_pos, max_pat_res,max_small_pat = 0,0,0,None
        for pat_res in np.arange(40,500,sizeres):#np.arange(10,500):
            #pat is 300 long, there are 300 pixels in vs
            #small_pat is pat_res long.
            small_pat = np.interp(np.linspace(0,29,int(pat_res)),np.linspace(0,29,20*6),pat)
            norm_small_pat=small_pat -np.mean(small_pat)
            norm_small_pat/=np.std(norm_small_pat)
            cross_cor = np.correlate(norm_vs,norm_small_pat)/len(small_pat)
            m = np.max(cross_cor)
            #print(m)
            if m>maxv:
                maxv = m
                max_pos = np.argmax(cross_cor)
                max_pat_res = pat_res
                max_small_pat = norm_small_pat
                max_cross_cor = cross_cor
        start = max_pos
        if max_small_pat is None:
            continue
        step = len(max_small_pat)/6
        items = np.linspace(start-step*2,start+step*9,12).astype(int)
        #print("MAX POS:")
        #print(max_pos)
        #plt.plot(norm_vs)
        #plt.plot(max_pos+np.arange(len(max_small_pat)),10*max_small_pat)
        #plt.vlines(items,0,10)
        raw = []
        for itstart,itend in zip(items[:-1]+int(step/4),items[1:]-int(step/4)):
            raw.append(np.mean(norm_vs[itstart:itend]))
        raw = np.array(raw)
        bindata = raw<np.mean(raw)
        
        if not np.all(bindata[2:7]==np.array([False,True,False,True,True])):
            #print("Search pattern error")
            continue
        #01 234567 89
        #DD 010110 DDP
        result = bindata[0]*8 + bindata[1]*4 + bindata[8]*2 + bindata[9]*1
        parity = bindata[10]*1
        print(result,np.min(np.abs(raw)))
        if np.min(np.abs(raw))<0.25:
            #print("raw value too small")
            continue
        if sum(bindata)%2==0:
            #print("Parity bit failed")
            continue
        #print(result)
        #np.set_printoptions(precision=1,suppress=True)
        #print(raw)
        #print(maxv)
        if optimum_result is not None:
            if maxv<optimum_result['maxv']:
                continue
        if np.abs(start-step*2)>20:
            continue
        #side info for optimum_result...
        dists = max_pos + np.arange(-2,10)*step
        xlocs = (np.cos(ang)*dists+markerloc[0]).astype(int)
        ylocs = (np.sin(ang)*dists+markerloc[1]).astype(int)
        optimum_result = {'code':result,'xlocs':xlocs,'ylocs':ylocs,'angle':ang,'maxv':maxv,'step':step,'start':start-step*2,'max_pat_res':max_pat_res,'max_pos':max_pos,'markerloc':markerloc}
    return optimum_result
    
