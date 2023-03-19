
def get_weighted_average(yields, v2s):

    weighted_v2 = 0
    weighted_yield = 0

    for i in range(len(yields)):
        weighted_v2 += yields[i] * v2s[i]
        weighted_yield += yields[i]

    return weighted_v2 / weighted_yield

assoc_pt_values = [1.5, 1.8, 2.4, 2.8, 3.2, 3.8]
trigger_pt_values = [4.5, 5.5]

lambda_v2_values = [0.063, 0.095, 0.12, 0.16, 0.155, 0.16]
associated_v2_values = [0.095, 0.11, 0.12, 0.125, 0.11, 0.105 ]
trigger_v2_values = [0.1, 0.075]

lambda_yield_values = [0.077, (0.055 + 0.040)/2, (0.020 + 0.014)/2, (0.0099 + 0.0059)/2, (0.0059 + 0.0026)/2, 0.00118]
associated_yield_values = [(0.49 + 0.37)/2, (0.23 + 0.18)/2, (0.057 + 0.046)/2, (0.025 + 0.021)/2, (0.013 + 0.011)/2, (0.008 + 0.007)/2]

trigger_yield_values = [0.001, 0.0005]



lambda_v2_central = get_weighted_average(lambda_yield_values[1:], lambda_v2_values[1:])
associated_v2_central = get_weighted_average(associated_yield_values[1:], associated_v2_values[1:])

lambda_v2_lowpt = get_weighted_average(lambda_yield_values[0:2], lambda_v2_values[0:2])
associated_v2_lowpt = get_weighted_average(associated_yield_values[0:2], associated_v2_values[0:2])

lambda_v2_highpt = get_weighted_average(lambda_yield_values[2:], lambda_v2_values[2:])
associated_v2_highpt = get_weighted_average(associated_yield_values[2:], associated_v2_values[2:])

trigger_v2 = get_weighted_average(trigger_yield_values, trigger_v2_values)


print("Lambda v2 central: ", lambda_v2_central)
print("Associated v2 central: ", associated_v2_central)
print("Lambda v2 lowpt: ", lambda_v2_lowpt)
print("Associated v2 lowpt: ", associated_v2_lowpt)
print("Lambda v2 highpt: ", lambda_v2_highpt)
print("Associated v2 highpt: ", associated_v2_highpt)





print("Trigger v2: ", trigger_v2)
