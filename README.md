# Heston-Model

The objective of the project is split into two parts. The first objective is to build a one
month implied price distribution for the S&P 500 market index, which is the
underlying asset in our case. Using this price distribution, it is our goal to derive an
independent "at-risk" benchmark implying about the downside tail risk exposure to
investors in the SPX market index and corroborate it with the counterpart value
obtained from the CBOE Skew Index. The second objective is to use the pricing
application of the Heston Stochastic Volatility model and run a Monte Carlo
simulation to price a 100%/130% 90-day up-and-out Call option of the S&P 500
Index.
