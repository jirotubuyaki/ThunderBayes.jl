var documenterSearchIndex = {"docs":
[{"location":"#ThunderBayes.jl","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"","category":"section"},{"location":"#News","page":"ThunderBayes.jl","title":"News","text":"","category":"section"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"Under development & test. Pleae wait the first alpha release.","category":"page"},{"location":"#Installation","page":"ThunderBayes.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"ThunderBayes.jl is available through GitHub (https://github.com/jirotubuyaki/ThunderBayes.jl/). If download from GitHub, you can use devtools by the commands:","category":"page"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"Pkg.add(url = \"https://github.com/jirotubuyaki/ThunderBayes.jl/\")","category":"page"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"Once the packages are installed, it needs to be made accessible to the current Julia session by the commands:","category":"page"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"> using ThunderBayes","category":"page"},{"location":"#Functions","page":"ThunderBayes.jl","title":"Functions","text":"","category":"section"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"Normal Bayesian Nonparametric Clustering","category":"page"},{"location":"#Plan","page":"ThunderBayes.jl","title":"Plan","text":"","category":"section"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"We aim to develop bayesian nonparametric package by Julia. Lots of relate scientific papers are pubilshed. We organize the papers and develop useful library. Firstly We will implement clustering methods and aim to release the first alpha version.","category":"page"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"The academic papers for this project are below: ","category":"page"},{"location":"#Clusterings","page":"ThunderBayes.jl","title":"Clusterings","text":"","category":"section"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"Nieto-Barajas, L. E., &amp; Contreras-Cristán, A. (2014). A bayesian nonparametric approach for time series clustering. Bayesian Analysis, 9(1). https://doi.org/10.1214/13-ba852 \nBeraha, M., Guglielmi, A., &amp; Quintana, F. A. (2021). The semi-hierarchical Dirichlet process and its application to clustering homogeneous distributions. Bayesian Analysis, 16(4). https://doi.org/10.1214/21-ba1278 \nHeller, K. A., &amp; Ghahramani, Z. (2005). Bayesian hierarchical clustering. Proceedings of the 22nd International Conference on Machine Learning  - ICML '05. https://doi.org/10.1145/1102351.1102389 \nBacallado, S., Favaro, S., Power, S., &amp; Trippa, L. (2021). Perfect sampling of the posterior in the hierarchical pitman–yor process. Bayesian Analysis, -1(-1). https://doi.org/10.1214/21-ba1269 \nPage, G. L., &amp; Quintana, F. A. (2015). Predictions based on the clustering of heterogeneous functions via shape and subject-specific covariates. Bayesian Analysis, 10(2). https://doi.org/10.1214/14-ba919 \nMüller Peter. (2015). Bayesian nonparametric data analysis. Springer. ","category":"page"},{"location":"#Inference","page":"ThunderBayes.jl","title":"Inference","text":"","category":"section"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"Ni, Y., Müller, P., Diesendruck, M., Williamson, S., Zhu, Y., &amp; Ji, Y. (2019). Scalable Bayesian nonparametric clustering and classification. Journal of Computational and Graphical Statistics, 29(1), 53–65. https://doi.org/10.1080/10618600.2019.1624366 \nCasella, G., Mengersen, K. L., Robert, C. P., &amp; Titterington, D. M. (2002). Perfect samplers for mixtures of distributions. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(4), 777–790. https://doi.org/10.1111/1467-9868.00360   \nNeal, R. M. (2003). Slice sampling. The Annals of Statistics, 31(3). https://doi.org/10.1214/aos/1056562461  ","category":"page"},{"location":"#Enviroment","page":"ThunderBayes.jl","title":"Enviroment","text":"","category":"section"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"Visual Studio Code - 1.68.1\njulia - 1.7.2","category":"page"},{"location":"#Contact","page":"ThunderBayes.jl","title":"Contact","text":"","category":"section"},{"location":"","page":"ThunderBayes.jl","title":"ThunderBayes.jl","text":"Please send suggestions at issues and report bugs to okadaalgorithm@gmail.com.","category":"page"},{"location":"Algorithms/crp_normal/#Bayesian-Nonparametric-Clustering-Algorithm","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"","category":"section"},{"location":"Algorithms/crp_normal/#Abstract","page":"Bayesian Nonparametric Clustering Algorithm","title":"Abstract","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Clustering is a scientific method which finds the clusters of data and many related methods are traditionally researched. Bayesian nonparametrics is statistics which can treat models having infinite parameters. Chinese restaurant process is used in order to compose Dirichlet process. The clustering which uses Chinese restaurant process does not need to decide the number of clusters in advance. This algorithm automatically adjusts it. Then, this package can calculate clusters in addition to entropy as the ambiguity of clusters.","category":"page"},{"location":"Algorithms/crp_normal/#Introduction","page":"Bayesian Nonparametric Clustering Algorithm","title":"Introduction","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Clustering is an analytical method in order to find the clusters of data and many related methods are proposed. K-means[1] and Hierarchical clustering[2] are famous algorithmic methods. Density-based clustering[3] is the method that finds clusters by calculating a concentration of data. In statistical methods, there are stochastic ways such as bayesian clustering[4]. However these methods need to decide the number of clusters in advance. Therefore if the data is both high dimensions and a complex, deciding the accurate number of clusters is difficult. Bayesian nonparametric method[5] composes infinite parameters by Dirichlet process[6]. Dirichlet process is the infinite dimensional discrete distribution that is composed by Stocastic processes like a chinese restaurant process (CRP)[7] or stick-breaking process[8]. CRP does not need to decide the number of clusters in advance. This algorithm automatically adjusts it. We implement the CRP Clustering and the method which calculates the entropy[9] into R package. Then, we explain the clustering model and how to use it in detail and execute simulation by example datasets.","category":"page"},{"location":"Algorithms/crp_normal/#Background","page":"Bayesian Nonparametric Clustering Algorithm","title":"Background","text":"","category":"section"},{"location":"Algorithms/crp_normal/#Chinese-Restaurant-Process","page":"Bayesian Nonparametric Clustering Algorithm","title":"Chinese Restaurant Process","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Chinese restaurant process is a metaphor looks like customers sit at a table in Chinese restaurant. All customers except for xi have already sat at finite tables. A new customer xi will sit at either a table which other customers have already sat at or a new table. A new customer tends to sit at a table which has the number of customers more than other tables. A probability equation is given by    ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"(Image: equa)","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"where n^i_k denotes the number of the customers at a table k except for i and α is a concentration parameter.","category":"page"},{"location":"Algorithms/crp_normal/#Markov-Chain-Monte-Carlo-Methods-for-Clustering","page":"Bayesian Nonparametric Clustering Algorithm","title":"Markov Chain Monte Carlo Methods for Clustering","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Markov chain Monte Carlo (MCMC) methods[10] are algorithmic methods to sample from posterior distributions. If conditional posterior distributions are given by models, it is the best way in order to acquire parameters from posterior distributions. The algorithm for this package is given by  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"i) ii) iterations continue on below:  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"i) Sampling z_i for each i (i = 1,2, ・・・,n)","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"(Image: equa)","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"where k is a k th cluster and i is a i th data. μ new and Σ new are calculated from all data.","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"ii) Calculating parameters for each k (k = 1,2, ・・・,∞)","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"(Image: equa)","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Iterations i) ii) continue by iteration number, and Σ k is a variance-covariance matrix of kth cluster. Cov is covariance. i and j are rows and columns’ number of Σ k ij and Σ new ij . First several durations of iterations which are called as burnin are error ranges. For that reason, burnin durations are abandoned.","category":"page"},{"location":"Algorithms/crp_normal/#Clusters-Entropy","page":"Bayesian Nonparametric Clustering Algorithm","title":"Clusters Entropy","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Entropy denotes the ambiguity of clustering. As a result of a simulation, data x i joins in a particular cluster. From the total numbers n k of the particular cluster k at the last iteration, a probability p k at each cluster k is calculated. The entropy equation is given by","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"(Image: equa)","category":"page"},{"location":"Algorithms/crp_normal/#API-Methods","page":"Bayesian Nonparametric Clustering Algorithm","title":"API Methods","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"ThunderBayes.ThunderBayes\ndata_check\ncrp_train\ncrp_predict\ncrp_visualize","category":"page"},{"location":"Algorithms/crp_normal/#ThunderBayes.ThunderBayes","page":"Bayesian Nonparametric Clustering Algorithm","title":"ThunderBayes.ThunderBayes","text":"ThunderBayes\n\nBayesian Nonparametric Clustering Algorithms\n\n\n\n\n\n","category":"module"},{"location":"Algorithms/crp_normal/#ThunderBayes.data_check","page":"Bayesian Nonparametric Clustering Algorithm","title":"ThunderBayes.data_check","text":"data_check(data)\n\nCheck whether data contain NaN value.\n\nExamples\n\n    ok = data_check(data)\n\nArguments\n\ndata::Array : the colums are the values of dimentions.\n\nReturn\n\nok::Bool : whether data contain NaN or not.\n\n\n\n\n\n","category":"function"},{"location":"Algorithms/crp_normal/#ThunderBayes.crp_train","page":"Bayesian Nonparametric Clustering Algorithm","title":"ThunderBayes.crp_train","text":"crp_train\n\nCalculate the clusters of data. MCMC Algorithms have been implemented. \n\nExamples\n\nburn_in = 5000\niteration = 50000\nmu = [0, 0, 0, 0, 0]\nsigma_table = Matrix{Float64}(1 * I, 5, 5)\nalpha = 1\nro_0 = 1\nautoparams = false\nresult, cluster_id = crp_train(data, burn_in, iteration, mu, sigma_table, alpha, ro_0, auto_params)\n\nArguments\n\ndata::Array : the colums are the values of dimentions.\nburn_in::Integer : an iteration integer of burn in.  burn in duration is abandoned.\niteration::Integer : an iteration integer.   \nmu::Array : a vector of center points of data. If data is 3 dimensions, a vector of 3 elements like \"[2, 4, 1]\".  \nsigma_table::Matrx : a numeric of table position variance.  \nalpha::Integer=1 : a numeric of a CRP concentration rate.  \nro_0::Integer=1 : a numeric of a CRP mu change rate.\nautoparams::Bool : automate sigma_table and mu for your data.\n\nReturn\n\nResult::Array : the mean values and variance-covariance matrix of the Clusters .\ncluster_id::Vector : the cluster id of data. \n\n\n\n\n\n","category":"function"},{"location":"Algorithms/crp_normal/#ThunderBayes.crp_predict","page":"Bayesian Nonparametric Clustering Algorithm","title":"ThunderBayes.crp_predict","text":"crp_predict(data, result)\n\nPrediction for new data.\n\nExamples\n\n    prediction_result = crp_preduct(data, cluster_id)\n\nArguments\n\ndata::vector : a data for prediction.\nresult::Array : a return variable of function crp_train() .\n\nReturn\n\nprediction_result::Array : the first column is cluster id and next colums are joined probabiilty of each clusters.\n\n\n\n\n\n","category":"function"},{"location":"Algorithms/crp_normal/#ThunderBayes.crp_visualize","page":"Bayesian Nonparametric Clustering Algorithm","title":"ThunderBayes.crp_visualize","text":"crp_visualize(data, cluster_id)\n\nVisualization of clustering results.\n\nExamples\n\n    crp_visualize(data, cluster_id)\n\nArguments\n\ndata::Array : a colums are the values of dimentions.\ncluster_id::Vector : a return variable of function crp_train() .\n\n\n\n\n\n","category":"function"},{"location":"Algorithms/crp_normal/#Simulation","page":"Bayesian Nonparametric Clustering Algorithm","title":"Simulation","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"We use dataset from Clustering basic benchmark( http://cs.joensuu.fi/sipu/datasets/ )[11]. If increase α parameter, new clusters tend to increase. burin_in iterations are abandoned. The result is plotted and each data joins in any cluster. The graph is given by below:","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"(Image: equa)","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Figure 1: Aggregation: Data is 788 elements and 2 dimentions. Parameters are set as alpha=1, burnin=100, iteration=1000.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"(Image: equa)  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Figure 2: 3 normal distribution: Data is 1000 elements and 2 dimentions. Parameters are set as alpha=0.5, burnin=100, iteration=1000.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"(Image: equa)  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Figure 3: 10 dimentional normal distributions: Data is generated from 10 dimentional normal distributions and parameters are set as alpha=1, burnin=100, iteration=1000.  ","category":"page"},{"location":"Algorithms/crp_normal/#Conclusions","page":"Bayesian Nonparametric Clustering Algorithm","title":"Conclusions","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"Chinese restaurant process clustering was implemented and explained how to use it. Computer resources are limited. Computer processing power is the most important problem. After this, several improvements are planed. Please send suggestions and report bugs to okadaalgorithm@gmail.com.","category":"page"},{"location":"Algorithms/crp_normal/#Acknowledgments","page":"Bayesian Nonparametric Clustering Algorithm","title":"Acknowledgments","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"This activity would not have been possible without the support of my family and friends. To my family, thank you for much encouragement for me and inspiring me to follow my dreams. I am especially grateful to my parents, who supported me all aspects.  ","category":"page"},{"location":"Algorithms/crp_normal/#References","page":"Bayesian Nonparametric Clustering Algorithm","title":"References","text":"","category":"section"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[1] Hartigan, J. A.; Wong, M. A. Algorithmas136: A k-means clustering algorithm . Journal of the Royal Statistical Society, Series C. 28 (1): 100–108. JSTOR 2346830, 1979.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[2] Rokach, Lior, and Oded Maimon. \"Clustering methods.\" Data mining and knowledge discovery handbook. Springer US, 2005. 321-352.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[3] Ester, Martin; Kriegel, Hans-Peter; Sander, Jörg; Xu, Xiaowei (1996). Simoudis, Evangelos; Han, Jiawei; Fayyad, Usama M. (eds.). A density-based algorithm for discovering clusters in large spatial databases with noise. Proceedings of the Second International Con ference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[4] John W Lau & Peter J Green (2007) Bayesian Model-Based Clustering Procedures, Journal of Computational and Graphical Statistics, 16:3, 526-558, DOI: 10.1198/106186007X238855.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[5] Muller Peter, et al. Bayesian Nonparametric Data Analysis. Springer, 2015.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[6] Ferguson, Thomas. Bayesian analysis of some nonparametric problems. Annals of Statistics. 1 (2): 209–230., 1973.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[7] Pitman, Jim. Exchangeable and partially exchangeable random partitions. Probability Theory and Related Fields 102 (2): 145–158., 1995.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[8] Broderick, Tamara, et al. “Beta Processes, Stick-Breaking and Power Laws.” Bayesian Analysis, vol. 7, no. 2, 2012, pp. 439–476., doi:10.1214/12-ba715.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[9] Elliott H. Lieb; Jakob Yngvason. The physics and mathematics of the second law of thermodynamics. Physics Reports Volume:310 Issue:1 1-96., 1999.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[10] Liu, Jun S. The collapsed gibbs sampler in bayesian computations with applications to a gene regulation problem. Journal of the American Statistical Association 89 (427): 958–966., 1994.  ","category":"page"},{"location":"Algorithms/crp_normal/","page":"Bayesian Nonparametric Clustering Algorithm","title":"Bayesian Nonparametric Clustering Algorithm","text":"[11] P. Fränti and S. Sieranoja K-means properties on six clustering benchmark datasets. Applied Intelligence, 48 (12), 4743-4759, December 2018.","category":"page"}]
}
