function d = b2i(x)
    if size(x,2)>64
        error('cannot currently translate bit vectors longer than 64 bits')
    end
    d=uint64(0);
    for i=1:size(x,2)
        d = bitset(d,i,x(1,i));
    end
end