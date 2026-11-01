/* Inspect published Markdown in an installed Chromium browser.
 * Install playwright-core outside the repository and expose it with NODE_PATH.
 * node tools/check_github_rendering.cjs --browser <exe> --output <external-dir>
 *      --ref <commit> --files <JSON array of repository-relative paths>
 * Screenshots are review artifacts, not figures or numerical evidence.
 */
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');
const {execFileSync} = require('node:child_process');
const {chromium} = require('playwright-core');
const root = path.resolve(__dirname, '..');
const args = Object.fromEntries(process.argv.slice(2).reduce((pairs, value, i, all) => {
  if (i % 2 === 0) pairs.push([value.replace(/^--/, ''), all[i + 1]]);
  return pairs;
}, []));

function sourceMathCounts(text) {
  let fence = null, display = 0, inline = 0;
  for (const line of text.split(/\r?\n/)) {
    const match = line.match(/^(`{3,}|~{3,})(.*)$/);
    if (match) {
      if (fence === null) {
        fence = match[1];
        if (match[2].trim() === 'math') display++;
      } else if (match[1][0] === fence[0] && match[1].length >= fence.length && !match[2].trim()) {
        fence = null;
      }
    } else if (fence === null) {
      inline += [...line.matchAll(/\$`[^`\n]+`\$/g)].length;
    }
  }
  return {display, inline};
}

async function main() {
  if (!args.browser || !args.output || !args.ref || !args.files) {
    throw new Error('Required: --browser --output --ref --files');
  }
  const output = path.resolve(args.output);
  if (!path.relative(root, output).startsWith('..') && !path.isAbsolute(path.relative(root, output))) {
    throw new Error('Keep browser profiles and screenshots outside the repository.');
  }
  fs.mkdirSync(output, {recursive: false});
  const files = JSON.parse(fs.readFileSync(args.files, 'utf8').replace(/^\uFEFF/, ''));
  const browser = await chromium.launch({executablePath: args.browser, headless: true});
  const context = await browser.newContext({viewport: {width: 1440, height: 1100}});
  const report = {
    checked_at: new Date().toISOString(), ref: args.ref,
    browser: await browser.version(), viewport: {width: 1440, height: 1100},
    playwright_version: require('playwright-core/package.json').version,
    scope: 'Published GitHub pages. DOM checks and screenshots; screenshots require visual review.',
    pages: []
  };
  try {
    for (const file of files) {
      const source = execFileSync('git', ['show', `${args.ref}:${file}`], {cwd: root});
      const text = source.toString('utf8');
      const expected = sourceMathCounts(text);
      const page = await context.newPage();
      const pageErrors = [];
      page.on('pageerror', e => pageErrors.push(e.message));
      let row;
      try {
        const url = `https://github.com/Griboedow/chaos-attractor-lab/blob/${args.ref}/${file}`;
        const response = await page.goto(url, {waitUntil: 'domcontentloaded', timeout: 60000});
        await page.locator('article.markdown-body').waitFor({timeout: 30000});
        await page.waitForFunction(() => [...document.querySelectorAll('math-renderer')].every(e =>
          e.querySelector('math,mjx-container') || /error|missing|undefined|unrecognized|not allowed/i.test(e.innerText)
        ), {}, {timeout: 30000}).catch(() => {});
        const body = page.locator('article.markdown-body');
        const data = await body.evaluate(a => {
          const maths = [...a.querySelectorAll('math-renderer')];
          return {
            heading_count: a.querySelectorAll('h1,h2,h3,h4').length,
            table_count: a.querySelectorAll('table:not(math table)').length,
            display_math_count: maths.filter(e => e.classList.contains('js-display-math')).length,
            total_math_count: maths.length,
            math_errors: maths.filter(e => !e.querySelector('math,mjx-container') || e.querySelector('merror,mjx-merror')).map(e => e.innerText),
            unsupported_labeled_rows: a.querySelectorAll('math mlabeledtr').length,
            images: [...a.querySelectorAll('img')].map(e => ({src: e.getAttribute('src'), loaded: e.complete && e.naturalWidth > 0})),
            width: a.clientWidth, scroll_width: a.scrollWidth,
            raw_display_delimiters: /(^|\n)\s*\$\$/.test(a.innerText),
            formula_symbols: maths.filter(e => e.classList.contains('js-display-math')).map(e => [...e.querySelectorAll('mo')].map(m => m.textContent).join(''))
          };
        });
        // Images can finish after math; wait for them before recording success.
        await page.waitForFunction(() => [...document.querySelectorAll('article.markdown-body img')].every(e => e.complete), {}, {timeout: 60000}).catch(() => {});
        data.images = await body.locator('img').evaluateAll(images => images.map(e => ({src: e.getAttribute('src'), loaded: e.complete && e.naturalWidth > 0})));
        const screenshot = file.replaceAll('/', '--').replace(/\.md$/, '') + '.png';
        await page.screenshot({path: path.join(output, screenshot), fullPage: true, timeout: 60000});
        row = {
          file, url, status: response.status(),
          source_sha256: crypto.createHash('sha256').update(source).digest('hex'),
          expected_display_math: expected.display, expected_inline_math: expected.inline,
          ...data, page_errors: pageErrors,
          screenshot, screenshot_sha256: crypto.createHash('sha256').update(fs.readFileSync(path.join(output, screenshot))).digest('hex')
        };
        row.passed = row.status === 200 && !row.math_errors.length && !row.raw_display_delimiters &&
          row.unsupported_labeled_rows === 0 &&
          row.display_math_count === expected.display &&
          row.total_math_count === expected.display + expected.inline && row.images.every(e => e.loaded);
      } catch (error) {
        row = {file, passed: false, error: String(error)};
      } finally {
        await page.close();
      }
      report.pages.push(row);
      fs.writeFileSync(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n');
      console.log(JSON.stringify({file, passed: row.passed, display: row.display_math_count,
        errors: row.math_errors || row.error, images: row.images?.length}));
    }
  } finally {
    await browser.close();
  }
  report.passed = report.pages.every(p => p.passed);
  fs.writeFileSync(path.join(output, 'report.json'), JSON.stringify(report, null, 2) + '\n');
  if (!report.passed) process.exitCode = 1;
}
main().catch(error => {console.error(error); process.exitCode = 1;});
