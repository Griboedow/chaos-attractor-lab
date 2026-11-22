export function makeAdvance(parameters) {
  const {
    pump, field_retention: retention, base_detuning: detuning, fast_shift: fast,
    memory_gain: gain, relaxation
  } = parameters;
  return function advance([x, y, memory], steps, firstObserved = steps + 1) {
    let low = Infinity, high = -Infinity;
    for (let step = 1; step <= steps; step++) {
      const intensity = x * x + y * y;
      const response = intensity / (1 + intensity);
      const phase = detuning + fast * response + memory;
      const c = retention * Math.cos(phase), s = retention * Math.sin(phase);
      const real = pump + x * c - y * s;
      y = x * s + y * c;
      x = real;
      memory = (1 - relaxation) * memory + relaxation * gain * response;
      if (step >= firstObserved) {
        low = Math.min(low, memory);
        high = Math.max(high, memory);
      }
    }
    return { state: [x, y, memory], low, high };
  };
}

export function donorMatches(donors, advance) {
  return donors.starts.flatMap((start, group) => {
    let state = advance(start, donors.burn).state;
    return donors.fields.filter(field => field.id.startsWith(group ? "U" : "L")).map((donor, i) => {
      state = advance(state, donors.stride).state;
      return donor.step === donors.burn + (i + 1) * donors.stride &&
             state.every((value, j) => value === donor.state[j]);
    });
  });
}
