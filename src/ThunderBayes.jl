# Bayesian Nonparametric Clustering Algorithm
# 2022/06/27 Masashi OKADA 
# okadaalgorithm@gmail.com
# MIT LICENSE

"""
    ThunderBayes
Bayesian Nonparametric Clustering Algorithms
"""
module ThunderBayes

using Distributions, Random
using Statistics, StatsBase
using LinearAlgebra
using DataFrames
using Plots,ColorSchemes

export data_check

"""
    data_check(data)
Check whether data contain NaN value.

# Examples
```julia-repl
    ok = data_check(data)
```

# Arguments

- `data::Array` : the colums are the values of dimentions.

# Return

- `ok::Bool` : whether data contain NaN or not.

"""
function data_check(data)
    data_length::UInt128 = size(data)[1]
    dim::UInt128 = size(data)[2]
    for i = 1 : data_length
        for j = 1 : dim
            if isnan(data[i, j])
                println("Error: Data contain NaN")
                return false
            end
        end
    end
    return true
end

export crp_train
"""
    crp_train
Calculate the clusters of data. MCMC Algorithms have been implemented. 

# Examples
```julia-repl
burn_in = 5000
iteration = 50000
mu = [0, 0, 0, 0, 0]
sigma_table = Matrix{Float64}(1 * I, 5, 5)
alpha = 1
ro_0 = 1
autoparams = false
result, cluster_result = crp_train(data, burn_in, iteration, mu, sigma_table, alpha, ro_0, auto_params)
```
# Arguments

- `data::Array` : the colums are the values of dimentions.
- `burn_in::Integer` : an iteration integer of burn in.  burn in duration is abandoned.
- `iteration::Integer` : an iteration integer.   
- `mu::Array` : a vector of center points of data. If data is 3 dimensions, a vector of 3 elements like "[2, 4, 1]".  
- `sigma_table::Matrx` : a numeric of table position variance.  
- `alpha::Integer=1` : a numeric of a CRP concentration rate.  
- `ro_0::Integer=1` : a numeric of a CRP mu change rate.
- `autoparams::Bool` : automate sigma_table and mu for your data.

# Return

- `result::Array` : the mean values and variance-covariance matrix of the Clusters .
- `cluster_result::Vector` : the cluster id of data. 

"""
function crp_train(data, burn_in, iteration, mu, sigma_table, alpha, ro_0,auto_params)

    rng = MersenneTwister(1234)
    Random.seed!(rng)

    data_length::UInt128 = size(data)[1]
    dim::UInt128 = size(data)[2]

    mu_0 = Array{Float64}(undef, 1, dim)
    sigma_new = Matrix{Float64}(1 * I, dim, dim)

    if auto_params
        mu_in = Array{Float64}(undef, dim)
        for i = 1 : dim
            mu_in[i,1] = mean(data[:, i])
            mean_i = mean(data[:, i])
            for j = i : dim
                mean_j = mean(data[:, j])
                data_sum = 0
                for n = 1 : data_length
                    data_sum = data_sum + (data[n, i] - mean_i)*(data[n, j] - mean_j)
                end
                sigma_new[i, j] = data_sum / data_length
                sigma_new[j, i] = sigma_new[i, j]
            end
        end
        mu_0 = mu_in
    else
        mu_0 = mu
        sigma_new = sigma_table      
    end
  
    data_k = Dict{Int64, Array{Float64}}()
    data_k_index = Dict{Int64, Dict{Int64, Int64}}()
    mu_k = Dict{Int64, Array{Float64}}()
    sigma_k = Array{Float64}(undef, 1, dim, dim)
    z = Matrix{UInt128}(undef, iteration, data_length)
    n_k = Array{Int64}(undef, 1)
    k_count = 1

    prob_k = Array{Float64}(undef, 2)
    prob_k_crp = Array{Float64}(undef, 2)

    add_m = Array{Float64}(undef, 1, dim)
    add_m = [0 * x for x in 1 : dim]
    add_v = Array{Float64}(undef, 1 , dim, dim)
    add_n = Array{Int64}(undef, 1)
    add_n = 1
    add_p = Array{Float64}(undef, 1)

    prob_min = 2.225074e-400 #min folat64 is 2.225074e-308
    prob_min_value = 1.0e-200

    println("Iteration: " , 1 , " Table_number: " , k_count)
    for i = 1 : data_length
        if i == 1
            z[1,1] = 1
            p =  MvNormal(mu_0, sigma_new)
            sample_in = rand(p, 1)
            mu_k[1] = sample_in[:, 1]
            n_k[1] = 1
            data_k[1] = data[i, :]
            data_k_index[1] = Dict(i => 1)
            sigma_k[1, :, :] = sigma_new[:, :]

        else
            prob_k_sum = 0
            for j = 1 : k_count
                prob_k[j] =  pdf(MvNormal(mu_k[j], sigma_new), data[i, :])
                if isnan(prob_k[j])
                    prob_k[j] = prob_min_value
                    println("Note 2 : prob_k[j] is NaN")
                end
                if prob_k[j] < prob_min
                    prob_k[j] = prob_min_value
                    println("Note 3 : prob_k[j] is Small")
                end
                prob_k_crp[j] = n_k[j] / (data_length + alpha)
                if isnan(prob_k_crp[j])
                    prob_k_crp[j] = prob_min_value
                    println("Note 5 : prob_k_crp[j] is NaN")
                end
                if prob_k_crp[j] < prob_min
                    prob_k_crp[j] = prob_min_value
                    println("Note 5 : prob_k_crp[j] is Small")
                end              
                prob_k[j] =  prob_k[j] * prob_k_crp[j]
                if isnan(prob_k[j])
                    prob_k[j] = prob_min_value
                    println("Note 8 : prob_k[j] is NaN")
                end
                if prob_k[j] < prob_min
                    prob_k[j] = prob_min_value
                    println("Note 9 : prob_k[j] is Small")
                end
                prob_k_sum =  prob_k_sum + prob_k[j]
            end
            p = MvNormal(mu_0, sigma_new)
            sample_in = rand(p, 1)
            mu_k[k_count + 1] = sample_in[:, 1]
            prob_k[k_count + 1] =  pdf(MvNormal(mu_k[k_count + 1], sigma_new), data[i, :])
            if isnan(prob_k[k_count + 1])
                prob_k[k_count + 1] = prob_min_value
                println("Note 11 : prob_k[k_count + 1] is NaN")
            end
            if prob_k[k_count + 1] < prob_min
                prob_k[k_count + 1] = prob_min_value
                println("Note 12 : prob_k[k_count + 1] is Small")
            end
            prob_k[k_count + 1] =  prob_k[k_count + 1] * (alpha / (data_length + alpha))

            if isnan(prob_k[k_count + 1])
                prob_k[k_count + 1] = prob_min_value
                println("Note 14 : prob_k[k_count + 1] is NaN")
            end
            if prob_k[k_count + 1] < prob_min
                prob_k[k_count + 1] = prob_min_value
                println("Note 15 : prob_k[k_count + 1] is Small")
            end
            prob_k_sum =  prob_k_sum + prob_k[k_count + 1]
            for j = 1 : k_count +  1
                prob_k[j] =  prob_k[j] / prob_k_sum     
            end
            if isnan(prob_k_sum)
                println("Note 16 prob_k_sum is NaN")
            end
            p = Categorical(prob_k)
            sample_in = rand(p, 1)
            sample = sample_in[1]
            if k_count + 1 == sample
                k_count = k_count + 1
                data_k[k_count] = data[i, :]
                mu_k[k_count] = add_m
                sigma_k = [sigma_k; add_v]
                n_k = [n_k; add_n]
                prob_k = [prob_k; add_p]
                prob_k_crp = [prob_k_crp; add_p]
                sigma_k[k_count, : , :] = sigma_new[:, :]
                mu_k[k_count] = data[i, :]
                z[1, i]  = sample
                data_k[k_count] = data[i, :]
                data_k_index[k_count] = Dict(i => 1)
            else
                z[1, i]  = sample
                n_k[sample] = n_k[sample] + 1
                data_k[sample] = vcat(data_k[sample], data[i, :])
                value_max = 0
                for (key_in, value_in) in data_k_index[sample]
                    if value_in > value_max
                        value_max = value_in
                    end
                end
                merge!(data_k_index[sample], Dict(i => value_max + 1))
            end
        end
        for j = 1 : k_count
            data_k_j = transpose(reshape(data_k[j],(dim,n_k[j])))
            if n_k[j] > dim
                mu_k_j = mu_k[j]
                flag_sigma_k = false
                for n = 1 : dim
                    mean_n = mean(data_k_j[:, n])
                    mu_k_j[n, 1] = mean_n 
                    for o = n : dim
                        mean_o = mean(data_k_j[:, o])
                        data_k_j_sum = 0
                        for p = 1 : n_k[j]
                            data_k_j_sum = data_k_j_sum + (data_k_j[p, n] - mean_n)*(data_k_j[p, o] - mean_o)
                        end
                        sigma_k[j, n, o] = (data_k_j_sum / n_k[j])
                        sigma_k[j, o, n] = sigma_k[j, n, o]
                        if sigma_k[j ,n, o] == 0
                            flag_sigma_k = true
                        end
                    end
                end
                mu_k[j] = mu_k_j
                if flag_sigma_k
                    sigma_k[j, : , :] = sigma_new[:, :]
                end
            else
                mu_k_j = mu_k[j]
                for n = 1 : dim
                    mean_n = mean(data_k_j[:, n])
                    mu_k_j[n, 1] = mean_n
                end
                mu_k[j] = mu_k_j
                sigma_k[j, :, :] = sigma_new[:, :]
            end
        end
    end
    for t = 2 : iteration
        for i = 1 : data_length
            key = data_k_index[z[t - 1, i]]
            data_k_tmp = copy(data_k[z[t - 1 ,i]])
            deleteat!(data_k_tmp, (((key[i] - 1) * 2 + 1) : ((key[i] - 1) * 2 + dim)))
            data_k_index_tmp = Dict{Int64, Int64}()
            for (key_in, value_in) in data_k_index[z[t - 1, i]]
                if key[i] < value_in 
                    value_in = value_in - 1
                elseif i == key_in
                    key_in = i
                    value_in = 0
                end
                merge!(data_k_index_tmp, Dict(key_in => value_in))
            end
            data_k[z[t - 1, i]] = data_k_tmp   
            data_k_index[z[t - 1, i]] = data_k_index_tmp  
            n_k[z[t - 1, i]] = n_k[z[t - 1, i]] - 1

            prob_k_sum = 0
            for j = 1 : k_count
                if n_k[j] == 0
                    prob_k[j] = 0
                else
                    prob_k[j] =  pdf(MvNormal(mu_k[j], sigma_k[j, :, :]), data[i, :])
                    if isnan(prob_k[j])
                        prob_k[j] = prob_min_value
                        println("Note 2 : prob_k[j] is NaN")
                    end
                    if prob_k[j] < prob_min
                        prob_k[j] = prob_min_value
                        println("Note 3 : prob_k[j] is Small")
                    end
                    prob_k_crp[j] = n_k[j] / (data_length + alpha)
                    if isnan(prob_k_crp[j])
                        prob_k_crp[j] = prob_min_value
                        println("Note 5 : prob_k_crp[j] is NaN")
                    end
                    if prob_k_crp[j] < prob_min
                        prob_k_crp[j] = prob_min_value
                        println("Note 6 : prob_k_crp[j] is Small")
                    end              
                    prob_k[j] =  prob_k[j] * prob_k_crp[j]

                    if isnan(prob_k[j])
                        prob_k[j] = prob_min_value
                        println("Note 8 : prob_k[j] is NaN")
                    end
                    if prob_k[j] < prob_min
                        prob_k[j] = prob_min_value
                        println("Note 9 : prob_k[j] is Small")
                    end
                    prob_k_sum =  prob_k_sum + prob_k[j]
                end
            end
            p = MvNormal(mu_0, sigma_new)
            sample_in = rand(p, 1)
            mu_k[k_count + 1] = sample_in[:, 1]
            prob_k[k_count + 1] =  pdf(MvNormal(mu_k[k_count + 1], sigma_new), data[i, :])
            if isnan(prob_k[k_count + 1])
                prob_k[k_count + 1] = prob_min_value
                println("Note 11 : prob_k[k_count + 1] is NaN")
            end
            if prob_k[k_count + 1] < prob_min
                prob_k[k_count + 1] = prob_min_value
                println("Note 12 : prob_k[k_count + 1] is Small")
            end
            prob_k[k_count + 1] =  prob_k[k_count + 1] * (alpha / (data_length + alpha))
            if isnan(prob_k[k_count + 1])
                prob_k[k_count + 1] = prob_min_value
                println("Note 14 : prob_k[k_count + 1] is NaN")
            end
            if prob_k[k_count + 1] < prob_min
                prob_k[k_count + 1] = prob_min_value
                println("Note 15 : prob_k[k_count + 1] is Small")
            end
            prob_k_sum =  prob_k_sum + prob_k[k_count + 1]
            for j = 1 : k_count +  1
                prob_k[j] =  prob_k[j] / prob_k_sum     
            end
            if isnan(prob_k_sum)
                println("Note 16 prob_k_sum is NaN")
            end
            p = Categorical(prob_k)
            sample_in = rand(p, 1)
            sample = sample_in[1] 
            if k_count + 1 == sample
                k_count = k_count + 1
                data_k[k_count] = data[i, :]
                mu_k[k_count] = add_m
                sigma_k = [sigma_k; add_v]
                n_k = [n_k; add_n]
                prob_k = [prob_k; add_p]
                prob_k_crp = [prob_k_crp; add_p]
                sigma_k[k_count, : , :] = sigma_new[:, :]
                mu_k[k_count] = data[i, :]
                z[t, i]  = sample
                data_k[k_count] = data[i, :]
                data_k_index[k_count] = Dict(i => 1)
            else
                z[t, i]  = sample
                n_k[sample] = n_k[sample] + 1
                data_k[sample] = vcat(data_k[sample], data[i, :])
                value_max = 0
                for (key_in, value_in) in data_k_index[sample]
                    if value_in > value_max
                        value_max = value_in
                    end
                end
                merge!(data_k_index[sample], Dict(i => value_max + 1))
            end
        end
        for j = 1 : k_count
            data_k_j = transpose(reshape(data_k[j],(dim,n_k[j])))
            if n_k[j] > dim
                mu_k_j = mu_k[j]
                flag_sigma_k = false
                for n = 1 : dim
                    mean_n = mean(data_k_j[:, n])
                    mu_k_j[n, 1] = mean_n
                    for o = n : dim
                        mean_o = mean(data_k_j[:, o])
                        data_k_j_sum = 0
                        for p = 1 : n_k[j]
                            data_k_j_sum = data_k_j_sum + (data_k_j[p, n] - mean_n)*(data_k_j[p, o] - mean_o)
                        end
                        sigma_k[j, n, o] = (data_k_j_sum / n_k[j])
                        sigma_k[j, o, n] = sigma_k[j, n, o] 
                        if sigma_k[j ,n, o] == 0
                            flag_sigma_k = true
                        end
                    end
                end
                mu_k[j] = mu_k_j
                if flag_sigma_k
                    sigma_k[j, :, :] = sigma_new[:, :]
                end
            else
                mu_k_j = mu_k[j]
                for n = 1 : dim
                    mean_n = mean(data_k_j[:, n])
                    mu_k_j[n, 1] = mean_n
                end
                mu_k[j] = mu_k_j
                sigma_k[j, :, :] = sigma_new[:, :]
            end
        end
        println("Iteration: " , t , " Table_number: " , k_count)
    end
    z_count = Array{Int64}(undef, data_length, k_count + 1)
    max = Array{Int64}(undef, data_length)
    max_value = Array{Int64}(undef, data_length)
    for i = 1 : data_length
        max[i] = 0
        max_value[i] = 0
        for j = 1 : k_count
            z_count[i, j] = 0
        end
        for t = burn_in : iteration
            z_count[i, z[t, i]] = z_count[i, z[t, i]] + 1
        end
        for j = 1 : k_count
            if n_k[j] != 0
                if max_value[i] <= z_count[i, j]
                    max[i] = j
                    max_value[i] = z_count[i, j]
                end
            end
        end
    end
    n_k_result = Array{Int64}(undef, k_count + 1)
    for j = 1 : k_count
        n_k_result[j] = 0
    end
    for i = 1 : data_length
        n_k_result[max[i]] = n_k_result[max[i]] + 1
    end

    count = 1
    for j = 1 : k_count
        if n_k_result[j] >= 1
            count = count + 1
        end
    end

    result_arry = Array{Float64}(undef, count - 1,  2 + dim + dim * dim)
    cluster_index = Dict{Int64 ,Int64}()
    count = 1
    count_cluter = 1
    for j = 1 : k_count
        if n_k_result[j] >= 1
            cluster_index = merge!(cluster_index, Dict(j => count_cluter))
            result_arry[count, 1] = count_cluter
            result_arry[count, 2] = n_k[j]
            mu_k_tmp = mu_k[j]
            count_cluter = count_cluter + 1
            for k = 1 : dim
                result_arry[count, 2 + k] = mu_k_tmp[k, 1]
            end
            count_in = 1
            for k = 1 : dim
                for l = 1 : dim
                    result_arry[count, 2 + dim + count_in] = sigma_k[j, k ,l]
                    count_in = count_in + 1
                end
            end
            count = count + 1
        end
    end

    cluster_result = Array{Int64}(undef, data_length)
    for i = 1 : data_length
        cluster_result[i] = cluster_index[max[i]]
    end

    for i = 1 : size(result_arry)[1]
        println("id: " , result_arry[i, 1])
        println("n_k: " , result_arry[i, 2] , " mu: " , result_arry[i, 3 : 2 + dim])

        println("sigma: " , result_arry[i, 2 + dim + 1 : 2 + dim + dim * dim])
    end
    println("cluster number: " , size(result_arry)[1])

    entropy_all = 0
    for i = 1 : size(result_arry)[1]
        entropy_all = entropy_all + (result_arry[i, 2] / data_length) * log2(result_arry[i, 2] / data_length)
    end
    println("entropy: " , -1 * entropy_all)
    
    return result_arry, cluster_result
end

export crp_visualize
"""
    crp_visualize(data, cluster_result)
Visualization of clustering results.

# Examples
```julia-repl
    crp_visualize(data, cluster_id)
```

# Arguments

- `data::Array` : a colums are the values of dimentions.
- `cluster_result::Vector` : a return variable of function crp_train() .

"""
function crp_visualize(data, cluster_result)
    data_dim = size(data)[2]
    gr()
    plt = []
    for i = 1 : data_dim
        for j = 1 : data_dim
            if i == data_dim
                x_label_str = string(j, " dim")
            else
                x_label_str = ""
            end
            if j == 1
                y_label_str = string(i , " dim")
            else
                y_label_str = ""
            end
            if i != j
                push!(plt, plot(data[:, i], data[:, j],
                    st =:scatter,
                    markersize = 2.4,
                    markerstrokewidth = 0.1,
                    markerstrokealpha = 0.0,
                    group = cluster_result,
                    color = cluster_result,
                    palette = :tab10,
                    xguidefontsize = 7,
                    yguidefontsize = 7,
                    xtickfontsize = 4,
                    ytickfontsize = 4,
                    xlabel = x_label_str,
                    ylabel = y_label_str,
                    legendfontsize = 4,
                    size = (1000 / data_dim, 1000 /data_dim),
                    legend = true))
            else   
                push!(plt, histogram(data[:, i],
                    bins = :scott,
                    palette = :tab10,
                    xguidefontsize = 7,
                    yguidefontsize = 7,
                    xtickfontsize = 4,
                    ytickfontsize = 4,
                    xlabel = x_label_str,
                    ylabel = y_label_str,
                    legend = false))          
            end
        end
    end
    plot(plt..., layout = (data_dim,data_dim), size = (1000, 1000))
end

export crp_predict
"""
    crp_predict(data, result)
Prediction for new data.

# Examples
```julia-repl
    prediction_result = crp_preduct(data, cluster_id)
```

# Arguments

- `data::vector` : a data for prediction.
- `result::Array` : a return variable of function crp_train() .


# Return

- `prediction_result::Array` : the first column is cluster id and next colums are joined probabiilty of each clusters.

"""
function crp_predict(data, result)

end
end
