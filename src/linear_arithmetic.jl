# addition, subtraction and scalar multiplication of symmetric functions

# addition of symmetric functions
function Base.:+(f::Sym,g::Sym)::Sym
    f_type=typeof(f)
    if f_type!=typeof(g)
        error("symmetric functions are not of the same type")
    end
    if f.degree!=g.degree
        error("symmetric functions are not of the same degree")   
    end
    p=min(length(f.coeffs),length(g.coeffs))
    return f_type(f.degree,f.coeffs[1:p]+g.coeffs[1:p])
end



# subtraction of symmetric functions
function Base.:-(f::Sym,g::Sym)::Sym
    f_type=typeof(f)
    f_type=typeof(f)
    if f_type!=typeof(g)
        error("symmetric functions are not of the same type")
    end
    if f.degree!=g.degree
        error("symmetric functions are not of the same degree")   
    end
    p=min(length(f.coeffs),length(g.coeffs))
    return f_type(f.degree,f.coeffs[1:p]-g.coeffs[1:p])
end




# multiplying symmetric functions with scalars
    
function Base.:*(f::Sym,c::Number)::Sym
    f_type=typeof(f)
    return f_type(f.degree,c.*f.coeffs)
end

function Base.:*(c::Number,f::Sym)::Sym
    f_type=typeof(f)
     return f_type(f.degree,c.*f.coeffs)
 end
    
