### ConstrainedGA

* Enter the command below to configure the anaconda environment.
<pre><code>conda create -c rdkit --name constrainedGA python=3.7</code></pre>

* Activate the anaconda environment.

<pre><code>conda activate constrainedGA</code></pre>

* Install all dependent packages with the command below.
<pre><code>conda install jupyter notebook
conda install -c conda-forge rdkit
conda install matplotlib
conda install seaborn
pip install selfies
</code></pre>

* paramenters(default)
<pre><code>--smiles_file 'data/zinc_dearom'.txt
--target_file 'data/logp_800.txt'
--seed 0
--population_size 100
--offspring_size 1000
--mutation_rate 0.5
--constraint_rate 0.1
--generations 20
--output_dir './output/constrained/'
--patience 5
--delta 0.6</code></pre>

* Example of Code Execution(delta = 0.6)
<pre><code>python constrained_opt.py</code></pre>

* Example of Code Execution(delta = 0.4)
<pre><code>python --delta 0.4 constrained_opt.py</code></pre>
