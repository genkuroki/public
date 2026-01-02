# The problem with p-values is how they’re used¹

* Andrew Gelman²
* 5 July 2013
* https://sites.stat.columbia.edu/gelman/research/unpublished/murtaugh2.pdf

1. Discussion of “In defense of P-values,” by Paul Murtaugh, for Ecology. We thank two reviewers for helpful comments and the National Science Foundation for partial support of this work.
2. Department of Statistics, Columbia University, New York, N.Y.

I agree with Murtaugh (and also with Greenland and Poole 2013, who make similar
points from a Bayesian perspective) that with simple inference for linear models, p-values
are mathematically equivalent to confidence intervals and other data reductions, there should
be no strong reason to prefer one method to another. In that sense, my problem is not with
p-values but in how they are used and interpreted.

Based on my own readings and experiences (not in ecology but in a range of social and
environmental sciences), I feel that p-values and hypothesis testing have led to much scientific
confusion by researchers treating non-significant results as zero and significant results as
real. In many settings I have found estimation rather than testing to be more direct. For
example, when modeling home radon levels (Lin et al. 1999), we constructed our inferences
by combining direct radon measurements with geographic and geological information. This
approach of modeling and estimation worked better than a series of hypothesis tests that
would, for example, reject the assumption that radon levels are independent of geologic
characteristics.

I have, on occasion, successfully used p-values and hypothesis testing in my own work,
and in other settings I have reported p-values (or, equivalently, confidence intervals) in
ways that I believe have done no harm, as a way to convey uncertainty about an estimate
(Gelman 2013). In many other cases, however, I believe that null hypothesis testing has led
to the publication of serious mistakes, perhaps most notoriously in the paper by Bem (2011),
who claimed evidence for extra-sensory perception (ESP) based on a series of statistically
significant results. The ESP example was widely recognized to indicate a crisis in psychology
research, not because of the substance of Bem’s implausible and unreplicated claims, but
because the research methods used to purportedly demonstrate the truth of these claims
were nothing but the standard null-hypothesis significance tests that are standard in so
many fields.

As researchers in medicine and psychology such as Ioannidis (2005), Simmons, Nelson,
and Simonsohn (2011), Yarkoni (2011), and Francis (2013) have discussed, the problem is not
merely with claims that are ridiculous on scientific grounds, but more broadly that many
statistically significant claims will be in error. Gelman and Weakliem (2009) discuss the
“statistical significance filter”: results that succeed in having a low p-value will inherently
yield overestimates of the magnitude of effects and comparisons (“Type M,” or magnitude,
errors) and are also likely to go in the wrong direction (“Type S,” or sign, errors).
The article under discussion reveals a perspective on statistics which, by focusing on
static data, is much different from mine. Murtaugh writes:

>Data analysis can be always be redone with different statistical tools. 
>The suitability of the data for answering a particular scientific question, however, cannot
>be improved upon once a study is completed. In my opinion, it would benefit 
>the science if more time and effort were spent on designing effective studies
>with adequate replication, and less on advocacy for particular tools to be used
>in summarizing the data.

I do not completely agree with this quotation, nor do I entirely agree with its implications.
First, the data in any scientific analysis are typically not set in stone, independent of the
statistical tools used in the analysis. Often I have found that the most important benefit
derived from a new statistical method is that it allows the inclusion of more data in drawing
scientific inferences. Here are some quick examples:

* Meta-analysis and hierarchical models allow partial pooling.
* Multinomial discrete-data regression models allow researchers to make fuller use of their measurements, going beyond the simple binary thresholding required for basic logistic regression.
* Multivariate methods such as factor analysis allow the use of multiple correlated mea- surements.
* Regularized regression methods such as lasso make it possible to include large numbers of predictors in regression models, much more than is possible using least squares methods for variable selection.

My second point of disagreement with the quotation above is in the implication that too
much time is spent on considering how to perform statistical inference. (Murtaugh writes
of “advocacy” but this seems to me to be a loaded term.) It is a well-accepted principle of
the planning of research that the design of data collection is best chosen with reference to
the analysis that will later be performed. We cannot always follow this guideline—once data
have been collected, they will ideally be made available for any number of analyses by later
researchers—but it still suggests that concerns of statistical methods are relevant to design.
Beyond this, as noted above, the choice of statistical method is not just about deciding how
to summarize “the data” but also influences what data are included in the analysis.

In conclusion, I share the long-term concern (see Krantz 1999, for a review) that the
use of p-values encourages and facilitates a sort of binary thinking in which effects and
comparisons are either treated as zero or are treated as real, and also an old-fashioned 
statistical perspective under which it is difficult to combine information from different sources.
The article under discussion makes a useful contribution by emphasizing that problems in
research behavior will not automatically be changed by changes in data reductions. The
mistakes that people make with p-values, could also be made using confidence intervals and
AIC comparisons, and I think it would be good for statistical practice to move forward from
the paradigm of yes/no decisions drawn from stand-alone experiments.

Hypothesis testing and p-values are so compelling in that they fit in so well with the
Popperian model in which science advances via refutation of hypotheses. For both theoretical
and practical reasons I am supportive of a (modified) Popperian philosophy of science in
which models are advanced and then refuted (Gelman and Shalizi 2013). But a necessary
part of falsificationism is that the models being rejected are worthy of consideration. If
a group of researchers in some scientific field develops an interesting scientific model with
predictive power, then I think it very appropriate to use this model for inference and to
check it rigorously, eventually abandoning it and replacing it with something better if it
fails to make accurate predictions in a definitive series of experiments. This is the form of
hypothesis testing and falsification that is valuable to me. In common practice, however,
the “null hypothesis” is a straw man that exists only to be rejected. In this case, I am
typically much more interested in the size of the effect, its persistence, and how it varies
across different situations. I would like to reserve hypothesis testing for the exploration of
serious hypotheses and not as in indirect form of statistical inference that typically has the
effect of reducing scientific explorations to yes/no conclusions.

## References

* Bem, D. J. 2011. Feeling the future: experimental evidence for anomalous retroactive influences on cognition and affect. Journal of Personality and Social Psychology 100:407–25.
* Francis, G. 2013. Replication, statistical consistency, and publication bias (with discussion). Journal of Mathematical Psychology.
* Gelman, A. 2013. P values and statistical practice. Epidemiology 24:69–72.
* Gelman, A., and Shalizi, C. 2013. Philosophy and the practice of Bayesian statistics (with discussion). British Journal of Mathematical and Statistical Psychology 66:8–18.
* Gelman A., and Weakliem, D. 2009. Of beauty, sex, and power: statistical challenges in estimating small effects. American Scientist 97:310316.
* Greenland S., and Poole, C. 2013. Living with P-values: resurrecting a Bayesian perspective on frequentist statistics. Epidemiology 24:62–68.
* Ioannidis, J. 2005. Why most published research findings are false. PLOS Medicine 2(8):e124.
* Krantz, D. H. 1999. The null hypothesis testing controversy in psychology. Journal of the American Statistical Association 44:1372–1381.
* Lin, C. Y., Gelman, A., Price, P. N., and Krantz, D. H. 1999. Analysis of local decisions using hierarchical modeling, applied to home radon measurement and remediation (with discussion). Statistical Science 14:305–337.
* Simmons J., Nelson L., and Simonsohn U. 2011. False-positive psychology: Undisclosed flexibility in data collection and analysis allow presenting anything as significant. Psy- chological Science 22:1359–1366.
* Yarkoni, T. (2011. The psychology of parapsychology, or why good researchers publishing good articles in good journals can still get it totally wrong. Citation Needed blog, 10 Jan. http://www.talyarkoni.org/blog/2011/01/10/the-psychology-of-parapsychology-or-why-good-researchers-publishing-good-articles-in-good-journals-can-still-get-it-totally-wrong/

--------------------

# p値の問題は、それがどのように使われているかにある¹

* **Andrew Gelman²**
* **2013年7月5日**
* [https://sites.stat.columbia.edu/gelman/research/unpublished/murtaugh2.pdf](https://sites.stat.columbia.edu/gelman/research/unpublished/murtaugh2.pdf)

1. Ecology誌に掲載された Paul Murtaugh による “In defense of P-values” へのディスカッション。本研究に対する有益なコメントをくださった2名の査読者に感謝する。また、本研究は米国国立科学財団（NSF）からの部分的支援を受けた。
2. コロンビア大学統計学部（ニューヨーク）

私は、単純な線形モデルにおける推測に関しては、p値が信頼区間やその他のデータ要約と数学的に同値であり、ある方法を他よりも強く好む理由はない、という点で Murtaugh（および同様の点をベイズ的観点から論じている Greenland と Poole 2013）に同意する。この意味において、私の問題は p値そのものにあるのではなく、それが **どのように使われ、どのように解釈されているか** にある。

私自身の読書経験や研究経験（生態学ではなく、社会科学や環境科学の広い分野において）に基づくと、p値や仮説検定は、「有意でない結果はゼロとして扱い、有意な結果は実在すると扱う」という研究者の姿勢を通じて、多くの科学的混乱を引き起こしてきたと感じている。多くの状況において、検定よりも推定のほうが、より直接的であると私は考えている。
例えば、住宅内ラドン濃度をモデル化した研究（Lin et al. 1999）では、直接測定されたラドン値と、地理学的・地質学的情報とを組み合わせて推論を構築した。このようなモデリングと推定のアプローチは、「ラドン濃度は地質学的特性と独立である」という仮定を棄却する、といった一連の仮説検定よりも、はるかにうまく機能した。

私は自分自身の研究において、p値や仮説検定をうまく使えたこともあるし、また別の状況では、推定値の不確実性を伝える手段として、p値（あるいは同値なものとしての信頼区間）を報告し、それが害を及ぼしたとは思っていない（Gelman 2013）。しかし他方で、帰無仮説検定が深刻な誤りの公表につながったと私が考える例も数多くある。その中でも特に悪名高いのが、Bem（2011）の論文であり、彼は一連の統計的に有意な結果に基づいて、超感覚的知覚（ESP）の証拠を主張した。
この ESP の例は、心理学研究における危機を示すものとして広く認識されたが、それは Bem の主張が非現実的で再現されていないからというよりも、むしろ、これらの主張の正しさを示すとされた研究方法が、多くの分野で標準的に用いられている、 **通常の帰無仮説有意性検定そのもの** であったからである。

Ioannidis（2005）、Simmons・Nelson・Simonsohn（2011）、Yarkoni（2011）、Francis（2013）といった医学・心理学分野の研究者が論じているように、問題は科学的に見て突飛な主張に限られない。より一般的に言えば、 **統計的に有意とされた主張の多くが誤っている** という点に問題がある。
Gelman と Weakliem（2009）は「統計的有意性フィルター」について論じている。すなわち、低い p値を達成した結果は、効果量や比較の大きさを体系的に過大評価する（Type M、すなわち magnitude error）傾向があり、さらに符号そのものが誤っている（Type S、すなわち sign error）可能性も高い。

本稿で議論されている論文は、統計学に対するある視点を示しているが、それは静的なデータに焦点を当てたものであり、私の立場とは大きく異なる。Murtaugh は次のように書いている。

> データ解析は、異なる統計的道具を用いて常にやり直すことができる。しかし、ある科学的問いに答えるうえでのデータの適切性は、ひとたび研究が完了してしまえば、それ以上改善することはできない。私の意見では、データを要約する際に用いる特定の道具を擁護することに時間を費やすよりも、十分な反復を備えた効果的な研究を設計することに、より多くの時間と労力が注がれるほうが、科学にとって有益であろう。

私はこの引用文にも、その含意にも、完全には同意しない。
第一に、科学的分析におけるデータは、通常、分析に用いられる統計的手法とは独立に、不変のものとして存在しているわけではない。新しい統計手法の最も重要な利点は、それによって **より多くのデータを推論に含めることが可能になる** 点にある、ということを私はしばしば経験してきた。簡単な例をいくつか挙げよう。

* メタ解析や階層モデルは、部分プーリングを可能にする。
* 多項離散データ回帰モデルは、単純なロジスティック回帰で必要とされる二値化を超えて、測定値をより完全に活用することを可能にする。
* 因子分析のような多変量手法は、相関した複数の測定を利用することを可能にする。
* lasso のような正則化回帰法は、変数選択に最小二乗法を用いる場合には不可能なほど多くの説明変数を、回帰モデルに含めることを可能にする。

第二の不一致点は、統計的推論の方法を検討することに「時間がかかりすぎている」という含意である（Murtaugh は「擁護（advocacy）」という言葉を使っているが、私にはこれは含みのある言い方に思える）。
研究計画において、データ収集の設計は、その後に行われる分析を念頭に置いて選ばれるべきである、という原則は広く受け入れられている。実際には、すでにデータが収集されており、それが後続の研究者によってさまざまな分析に用いられる、という状況も少なくないが、それでもなお、この原則は、統計的方法に関する配慮が研究設計にとって重要であることを示唆している。
さらに言えば、前述したように、統計的方法の選択は単に「データ」をどう要約するかの問題にとどまらず、 **どのデータが分析に含められるか** そのものに影響を与える。

結論として、私は（レビューとしては Krantz 1999 を参照）p値の使用が、効果や比較を「ゼロか、実在か」という二分的思考を助長し、また、異なる情報源からの情報を統合することが難しい、旧式の統計的視点を促進している、という長年の懸念を共有している。
本稿で議論されている論文は、研究行動の問題は、データ要約の方法を変えるだけでは自動的には解決しない、という点を強調しており、その点で有益な貢献をしている。p値で犯される誤りは、信頼区間や AIC 比較を用いても犯されうるものである。私は、統計的実践が、単発の実験から yes/no の決定を引き出すというパラダイムを超えて前進することが望ましいと考えている。

仮説検定と p値は、「科学は仮説の反証によって進歩する」というポパー的モデルに非常によく適合するという意味で、きわめて魅力的である。理論的・実践的理由の両面から、私は、モデルが提案され、それが反証されるという（修正された）ポパー主義的科学観を支持している（Gelman and Shalizi 2013）。
しかし、反証主義にとって不可欠なのは、棄却されるモデルが **検討に値するもの** であるという点である。ある科学分野の研究者集団が、予測力を持つ興味深い科学モデルを構築したのであれば、そのモデルを推論に用い、厳密に検証し、決定的な一連の実験で正確な予測ができないことが示されたならば、それを放棄して、より良いものに置き換えるのは、きわめて適切である。
これこそが、私にとって価値ある仮説検定と反証の姿である。しかし、一般的な実践においては、「帰無仮説」は、棄却されるためだけに存在する藁人形にすぎないことが多い。この場合、私が関心を持つのは、効果の大きさ、その持続性、そして状況によってそれがどのように変化するかである。
私は、仮説検定を、真剣に検討すべき仮説の探究のために用いるべきであり、科学的探究を yes/no の結論に還元してしまう、間接的な統計的推論の形式として用いるべきではないと考えている。

## 参考文献

* Bem, D. J. 2011. Feeling the future: experimental evidence for anomalous retroactive influences on cognition and affect. *Journal of Personality and Social Psychology* 100:407–425.
* Francis, G. 2013. Replication, statistical consistency, and publication bias (with discussion). *Journal of Mathematical Psychology*.
* Gelman, A. 2013. P values and statistical practice. *Epidemiology* 24:69–72.
* Gelman, A., and Shalizi, C. 2013. Philosophy and the practice of Bayesian statistics (with discussion). *British Journal of Mathematical and Statistical Psychology* 66:8–18.
* Gelman, A., and Weakliem, D. 2009. Of beauty, sex, and power: statistical challenges in estimating small effects. *American Scientist* 97:310–316.
* Greenland, S., and Poole, C. 2013. Living with P-values: resurrecting a Bayesian perspective on frequentist statistics. *Epidemiology* 24:62–68.
* Ioannidis, J. 2005. Why most published research findings are false. *PLOS Medicine* 2(8):e124.
* Krantz, D. H. 1999. The null hypothesis testing controversy in psychology. *Journal of the American Statistical Association* 44:1372–1381.
* Lin, C. Y., Gelman, A., Price, P. N., and Krantz, D. H. 1999. Analysis of local decisions using hierarchical modeling, applied to home radon measurement and remediation (with discussion). *Statistical Science* 14:305–337.
* Simmons, J., Nelson, L., and Simonsohn, U. 2011. False-positive psychology: undisclosed flexibility in data collection and analysis allow presenting anything as significant. *Psychological Science* 22:1359–1366.
* Yarkoni, T. 2011. The psychology of parapsychology, or why good researchers publishing good articles in good journals can still get it totally wrong. *Citation Needed* blog, 2011年1月10日。 [link](http://www.talyarkoni.org/blog/2011/01/10/the-psychology-of-parapsychology-or-why-good-researchers-publishing-good-articles-in-good-journals-can-still-get-it-totally-wrong/)

