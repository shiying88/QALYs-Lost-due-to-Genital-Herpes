import numpy as np
import SimPy.RandomVariateGenerators as RVG


class RandomValueGenerator:
    """ parent class """

    def __init__(self, note):
        self.note = note  # explanation of the variable


class BetaValueGenerator(RandomValueGenerator):
    """ variables follow beta distribution """

    def __init__(self, mean, stdv, note):
        super().__init__(note)

        fit_output = RVG.Beta.fit_mm(mean=mean, st_dev=stdv)
        self.rvg = RVG.Beta(a=fit_output['a'], b=fit_output['b'])

    def sample(self, rng):
        """ draw a random sample from parameterized beta distribution """
        return self.rvg.sample(rng)


class LogNormalValueGenerator(RandomValueGenerator):
    """ variables follow log-normal distribution """

    def __init__(self, note, mean=None, stdv=None, log_mean=None, log_std=None):
        super().__init__(note)

        # x follows a LogNormal distribution -> log(x) follows a normal distribution.
        # fit_mm(mean, st_dev), mean and st_dev are attributes for the observed sample (i.e., x)
        if mean is not None:
            fit_output = RVG.LogNormal.fit_mm(mean=mean, st_dev=stdv)
            self.rvg = RVG.LogNormal(mu=fit_output['mu'], sigma=fit_output['sigma'])
        # fit_output['mu'] and fit_output['sigma'] are mean and StDev for the underlying normal distribution (log(x))
        elif log_mean is not None:
            self.rvg = RVG.LogNormal(mu=log_mean, sigma=log_std)
        else:
            raise ValueError('wrong input of mean and st dev for log-normal distribution')

    def sample(self, rng):
        """ draw a random sample from parameterized log-normal distribution """
        return self.rvg.sample(rng)


class GammaValueGenerator(RandomValueGenerator):
    """ variables follow gamma distribution """

    def __init__(self, mean, stdv, note):
        super().__init__(note)

        fit_output = RVG.Gamma.fit_mm(mean=mean, st_dev=stdv)
        self.rvg = RVG.Gamma(a=fit_output['a'], scale=fit_output['scale'])

    def sample(self, rng):
        """ draw a random sample from parameterized gamma distribution """
        return self.rvg.sample(rng)


class DirichletValueGenerator(RandomValueGenerator):
    """ variables follow dirichlet distribution """

    def __init__(self, N, prob_list, note):
        super().__init__(note)
        self.N = N  # population size
        self.prob_list = prob_list  # list of probabilities for every event
        self.note = note  # explanation
        self.events_num = []
        for i in range(len(self.prob_list)):
            event_num = self.N * self.prob_list[i]
            self.events_num.append(event_num)

    def sample(self, rng):
        """ draw a random sample from parameterized dirichlet distribution """
        # return np.random.dirichlet(alpha=self.events_num)
        return rng.dirichlet(alpha=self.events_num)
