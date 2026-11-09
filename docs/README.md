# oarfish documentation

The source of the oarfish website, https://combine-lab.github.io/oarfish/, built with
[Astro](https://astro.build) and its [Starlight](https://starlight.astro.build) documentation theme.

Pages are Markdown (or MDX) files under `src/content/docs/`; the sidebar is set in
`astro.config.mjs`. Links between pages include the site's `/oarfish` base, e.g.
`[Output files](/oarfish/reference/output/)`.

To preview locally (Node.js 22.23.4 or later):

```sh
cd docs
npm ci
npm run dev       # live preview at http://localhost:4321/oarfish/
npm run build     # the static site, in dist/
```

Pushes to `main` or `dev` that touch `docs/` build the site and publish it to the `gh-pages`
branch (`.github/workflows/deploy-docs.yml`); pull requests that touch `docs/` build it only.
