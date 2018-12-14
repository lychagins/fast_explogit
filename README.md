# fast_explogit
Fast mixed exploded logit with agent-specific choice sets.

## Usage ##
`[logl, grad] = explogit(beta, X, nskipped, nlisted)`

or

`[logl, grad] = explogit(beta, X, nskipped, nlisted, weight)`

### Arguments ###
`beta`: Vector of logit coefficients. Type: **double**.

`X`: Matrix of covariates. Columns correspond to covariates, rows correspond to agent-choice pairs. Type: **double**.

`nskipped`: Number of skipped choices for each agent. Type: **uint16**.

`nlisted`: Number of listed choices for each agent.
The number of listed choices, `nlisted`, and the total size of the choice set, `nlisted + nskipped`, are assumed to be exogenous. Type: **uint16**.

`weight`: _Optional argument_. Agent `i` in the sample represents `weight(i)` agents in the population. Default: the vector of ones. Type: **double**.

### Return values ###
`logl`: Loglikelihood function. Type: **double**.

`grad`: Gradient of `logl`. Type: **double**.

### Example ###

Suppose that agent A chooses 2 most preferred items from a choice set {1, 3, 4, 6}. Agent B chooses 1 most preferred item from {2, 5, 6}. A's first and second best choices are 3 and 6 respectively, while B's first best is 2.

A's and B's choice sets and stated preferences determine the ordering of rows in `X`. Rows 1–4 represent the choice set of A, in the reverse preference order. That is, covariates of the most preferred item, 3, are placed to row 4; the second-best choice, 6, goes to row 3. Rows 1 and 2 correspond to choices 1 and 4 (their relative order is unimportant). The block in rows 5–7 correspond to choices of agent B: item 2 is in row 7, items 5 and 6 are in rows 5–6.