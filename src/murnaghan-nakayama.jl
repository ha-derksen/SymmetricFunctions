# murnaghan-nakayama rule


# murnaghan_nakayama_rule1(lambda,k) computes the Murnaghan Nakayama rule for s_(lambda)*p_(k).
# The output is an integer vector (Int32). 
# if q appears in the vector, then the |q|-th partition of n+k appears in s_(lambda)*p_(k)
# it appears with a + sign if q>0 and with a - sign if q<0
function murnaghan_nakayama_rule1(lambda::Vector{<:Integer},k::Integer)::Vector{<:Integer}
    max_MN1_output=sum(lambda)+2*k-1
    lambda2=[lambda;zeros(eltype(lambda),k)]
    i=1 
    t=0
    MN1_output=zeros(Int,max_MN1_output)
    for j in eachindex(lambda2)
        while k+lambda2[j]+i<=lambda2[i]+j 
            i+=1
        end
        if  k+lambda2[j]+i>lambda2[i]+j  && (i==1 || k+lambda2[j]+i<=lambda2[i-1]+j)
                mu=copy(lambda2)       
                mu[i]=lambda2[j]+k+i-j
                mu[i+1:j]=lambda2[i:j-1].+1  
                t+=1
                MN1_output[t]=(-1)^(i-j)*partition2rank(mu)
        end
    end
    return MN1_output[1:t]
end



    



# murnaghan_nakayama_rule1(n::integer,i::Integer,k::Integer) computes murnaghan_nakayama_rule1(lambda,k) where k is the i-th partition of n.
function murnaghan_nakayama_rule1(n::Integer,i::Integer,k::Integer)::Vector{<:Integer}
    if k==0
        return [Int(i)]
    end
    return murnaghan_nakayama_rule1(rank2partition(n,i),k)
end







# 'murnaghan_nakayama_rule2(lambda,k)' computes the Murnaghan Nakayama rule when applying differentiation with respect 
# to p_(k) to s_(lambda). The output is an integer vector. 
# If the integer q appears in the vector, then the |q|-th partition of n+k
# appears in the product s_(lambda)*p_(k). It appears with a + sign if q>0 and with a - sign if q<0.
function murnaghan_nakayama_rule2(lambda::Vector{<:Integer},k::Integer)::Vector{<:Integer}
    max_MN2_output=Int(floor((sqrt(k^2+8*sum(lambda))-k)/2+.001))
    n=length(lambda)
    j=1
    t=0
    MN2_output=zeros(Int,max_MN2_output)
    for i in 1:n
        while (j<n) && (j<i || lambda[i]-k+j-i+1<=lambda[j+1])
            j+=1
        end
        if (lambda[i]-k+j-i+1<=lambda[j]) && lambda[i]-k+j-i+1>0
            mu=copy(lambda)
            mu[i:j-1]=lambda[i+1:j].-1
            mu[j]=lambda[i]-k+j-i
#            mu=[lambda[1:i-1];lambda[i+1:j].-1;lambda[i]-k+j-i;lambda[j+1:n]]
            t+=1
            MN2_output[t]=(-1)^(i-j)*partition2rank(mu)
#            append!(MN2_output,(-1)^(i-j)*partition2rank(mu))
        end
    end
    return MN2_output[1:t]
end




# murnaghan_nakayama_rule2(n::integer,i::Integer,k::Integer) computes murnaghan_nakayama_rule2(lambda,k) where k is the i-th partition of n.
function murnaghan_nakayama_rule2(n::Integer,i::Integer,k::Integer)::Vector{<:Integer}
    if k==0
        return [Integer(i)]
    end
    return murnaghan_nakayama_rule2(rank2partition(n,i),k)
end


