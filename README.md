# CNNClassifier
Convolutional Neural Net for Sequence Classification

In this project, we are interested in finding the motif for TOP2B. It has been discovered that CTCF and RAD21 binding sites overlapped, and CTCF+RAD21 sites are biologically important. Recently it has been shown that TOP2B also associates with CTCF and RAD21. 

The problem is transformed into a binary classification problem, where the two classes are:
* double sites: CTCF + RAD21
* triple sites: CTCF + RAD21 + TOP2B 

## Models:
* machining learning methods: logistic regression, random forest, boosting, SVM (Support Vector Machine)
* CNN 

## Content:
* classify double and triple sites (ML & CNN)
* find the motif for TOP2B sites through model visualization and interpretation
* study the association between triple sites and mutations
