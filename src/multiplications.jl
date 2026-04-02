# operations & (S_n tensor product) and * (GL_n rep. product)



# f&g computes the (tensor) product of f and g viewed as characters of the symmetric group S_n.
# f and g must be the same type of symmetric function and should have the same degree n.
function Base.:&(f::Sym,g::Sym)::Sym
    f_type=typeof(f)
    if f_type!=typeof(g)
        error("symmetric functions are not of the same type")
    end
    nf=f.degree
    ng=g.degree
    if nf!=ng
        error("symmetric functions are not of the same degree")
    end
    f_length=length(f.coeffs)
    g_length=length(g.coeffs)
    min_length=min(f_length,g_length)
    fc=SymC(f)
    gc=SymC(g)
    h=f_type(SymC(nf,fc.coeffs[1:min_length].*gc.coeffs[1:min_length]))
    h_coeff_type=common_coeff_type(f,g)
    return f_type(nf,Vector{h_coeff_type}(h.coeffs))
end


# f*g computes the tensor product decomposition of two representations of GL_n (where n is arbitrarily large)
function Base.:*(f::Sym,g::Sym)::Sym
    type_of_f=typeof(f)
    if type_of_f!=typeof(g)
          error("symmetric functions are not of the same type")
    end
    fp=SymP(f)
    gp=SymP(g)
    fdeg=f.degree
    gdeg=g.degree
    n=fdeg+gdeg
    CoeffType_h=common_coeff_type(f,g)
    CoeffType_h_p=common_coeff_type(fp,gp)
    h_p=SymP(n,zeros(CoeffType_h_p,n_partitions(n)))
    for i in eachindex(fp.coeffs)
        if fp.coeffs[i]!=0
            for j in eachindex(gp.coeffs)
                if gp.coeffs[j]!=0
                    r=partition2rank(merge_partitions(rank2partition(fdeg,i),rank2partition(gdeg,j)))
                    h_p.coeffs[r]+=fp.coeffs[i].*gp.coeffs[j]
                end
            end
        end
    end
    h=type_of_f(h_p)
    return type_of_f(n,Vector{CoeffType_h}(h.coeffs))
end

